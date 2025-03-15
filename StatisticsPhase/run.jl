import Pkg
Pkg.activate(Base.current_project())
Pkg.instantiate()

# load BLAS / MKL
use_MKL = parse(Bool, ARGS[1])
if use_MKL
    using MKL
    println("Using MKL")
    lib_string="_MKL"
else
    println("Using OPENBLAS")
    lib_string="_BLAS"
end

using StatsBase, Random, FileIO, LinearAlgebra, Hwloc, MLUtils
using Base.Threads: nthreads, @threads, @spawn
using TimeBinEncoding

Julia_threads = nthreads()
BLAS.set_num_threads(1)
BLAS_threads = BLAS.get_num_threads()
println("Found $(num_physical_cores()) physical and $(num_virtual_cores()) logical cores.")
println("Starting $Julia_threads Julia Threads and $BLAS_threads BLAS Threads.")

const n_resamples = parse(Int, ARGS[2]) # Number of repetitions of a full simulated experiment

const N = parse(Int, ARGS[3]) # Number of time bins
const N_samples = parse(Int, ARGS[4]) # Total number of samples per experiment
const ϵ_angles = parse(Float64, ARGS[5]) * π # ±ϵ error in the angles of the central coupler

const ϵ = parse(Float64, ARGS[6]) # dephasing strength
# 0.06905 for N = 2 | 0.054836 for N = 4 | 0.05215 for N = 8 | 0.05152 for N = 16

@assert ARGS[7][1] == '[' && ARGS[7][end] == ']'

const n_samples_phase_estimation_each_arr = parse.(Int, split(chop(ARGS[7], head=1),','))
    # Number of samples used for each phase-estimation DTQW

Random.seed!(8675309) # / Jenny dont change your number

function state_prep(
    ρ, angles_phase_real, angles_phase_imag, n_samples_pop, n_samples_phase_estimation_each
)
    # prepare the phase-corrected intitial state and compute all random/noisy quantities
    pops_init_sampled = populations(ρ, n_samples_pop)
        # sample initial-state time-bin populations
    pops_fs_real, pops_fs_imag = fs_pop_phase_estimation(
        ρ, angles_phase_real, angles_phase_imag
    ) # calculate the noisy final-state populations of the two phase estimation DTQWs
    pops_fs_real_sampled = sample_populations(
        pops_fs_real, n_samples_phase_estimation_each; unity=false
    ) # sample the final-state populations of the real-part phase-estimation DTQW
    pops_fs_imag_sampled = sample_populations(
        pops_fs_imag, n_samples_phase_estimation_each; unity=false
    ) # sample the final-state populations of the imaginary-part phase-esimation DTQW
    relative_phases_auto = initial_state_phase_estimation(
        pops_init_sampled, pops_fs_real_sampled, pops_fs_imag_sampled
    ) # compute the relative phases of the initial state
    ρ_corrected = phase_on_density_matrix(ρ, -1 * relative_phases_auto)
        # correct the initial state with the inverse phases
    return pops_init_sampled, pops_fs_real, pops_fs_imag, pops_fs_real_sampled, pops_fs_imag_sampled, ρ_corrected
end

function eval_single(
    N, ρ_corrected, pops_init_sampled, angles_single, angles_single_noisy, n_samples_single_setting
)
    # evaluate the fidelity using the single-setting DTQW scheme
    j_out_single = j_out_single_setup(N) # compute final-state projector indices
    j_contr_tuples =  correlated_short_bins_tuples(N; extract_diagonal = true)
        # indice tuples for contributing coherences, including the diagonal terms
    pop_fs_coh_single = explicit_fs_pop(
        ρ_corrected, j_out_single, angles_single_noisy
    ) # calculate the noisy final-state populations of the single-setting DTQW
    pop_fs_coh_single_sampled = sample_populations(
        pop_fs_coh_single, n_samples_single_setting
    ) # sample the final-state populations of the single-setting DTQW
    fidelity = coherence_extraction(
        N, j_out_single, pops_init_sampled, pop_fs_coh_single_sampled, angles_single, j_contr_tuples
    ) # compute the fidelity

    return pop_fs_coh_single_sampled, fidelity
end

Ψ_mes = insert_initial_state(correlated_timebin_state(fill(1 / sqrt(N), N)))
    # maximally entangled state wave function
angles_single = angles_single_setup(N)
    # angles for the single-setting DTQW
ρ_perfect_phase_correction = density_matrix_dephased(Ψ_mes, ϵ)
    # dephased density matrix, but with no relative phases, acting as a benchmark for the phase correction

for (idx_sample, n_sample_phase) in enumerate(n_samples_phase_estimation_each_arr)
    # iterate over the number of samples used for each phase-estimation DTQW

    n_samples_phase_estimation_each = n_sample_phase

    n_samples_populations = N_samples / 4 # number of samples used for the initial-state populations
    n_samples_compound_each =
        (N_samples - n_samples_populations - 2 * n_samples_phase_estimation_each) / (N - 1)
        # number of samples used for each of the N-1 compound-setting DTQWs

    n_samples_single_setting = (N - 1) * n_samples_compound_each
        # number of samples used for the one single-setting DTQW

    tasks_per_thread = 1
    # customize this as needed. More tasks have more overhead, but better load balancing

    n_threads_resampling = Int(Base.Threads.nthreads())
    data_chunks = chunk(1:n_resamples, tasks_per_thread * n_threads_resampling)
    # divides in potentially less than n_threads_resampling chunks, if no advantage is gained from more tasks
    tasks = map(data_chunks) do package
        # divide all resamples into chunks and distribute them to the threads
        Base.Threads.@spawn begin
            chunk_size = length(package)
            local fidelity_sampled_single_deph = zeros(chunk_size)
            local fidelity_sampled_compound_deph = zeros(chunk_size)
            local mes_fidelity_deph = zeros(chunk_size)

            local fidelity_sampled_single_pure = zeros(chunk_size)
            local fidelity_sampled_compound_pure = zeros(chunk_size)
            local mes_fidelity_pure = zeros(chunk_size)

            local phase_measurement_pure = zeros(N - 1, 4, chunk_size) # pop_real, pop_imag, pop_real_sampled, pop_imag_sampled
            local phase_measurement_deph = zeros(N - 1, 4, chunk_size) # 	-||-

            local pop_init_sampled_pure_arr = zeros(4 * N^2, chunk_size)
            local pop_init_sampled_deph_arr = zero(pop_init_sampled_pure_arr)

            local pop_fs_coh_noisy_single_sampled_deph_arr = zeros(chunk_size)
            local pop_fs_coh_noisy_single_sampled_pure_arr = zero(pop_fs_coh_noisy_single_sampled_deph_arr)

            local pop_fs_coh_noisy_compound_sampled_deph_arr = zeros(N - 1, chunk_size)
            local pop_fs_coh_noisy_compound_sampled_pure_arr = zero(pop_fs_coh_noisy_compound_sampled_deph_arr)

            local angles_noisy_compound_arr = []
            local angles_noisy_single_arr = []
            local angles_noisy_phase_real_arr = []
            local angles_noisy_phase_imag_arr = []

            for i in 1:chunk_size
                # iterate over each independent experimental run
                wf_coeffs = cis.(2 * rand(N) * π)
                    # random phase coefficients for the initial state
                Ψ_init = insert_initial_state(correlated_timebin_state(wf_coeffs))
                ρ_pure = density_matrix(Ψ_init)
                    # pure random-phase intitial state
                ρ_mixed = density_matrix_dephased(Ψ_init, ϵ)
                    # dephased random-phase initial state


                angles_noisy_phase_real = angles_phase_estimation(N, ϵ_angles)
                angles_noisy_phase_imag = angles_phase_estimation(N, ϵ_angles)
                    # noisy angles for the phase estimation DTQWs

                pops_init_sampled_pure, pops_fs_real_pure, pops_fs_imag_pure, pops_fs_real_sampled_pure, pops_fs_imag_sampled_pure, ρ_corrected_pure =
                    state_prep(ρ_pure, angles_noisy_phase_real, angles_noisy_phase_imag, n_samples_populations, n_samples_phase_estimation_each)
                pops_init_sampled_deph, pops_fs_real_deph, pops_fs_imag_deph, pops_fs_real_sampled_deph, pops_fs_imag_sampled_deph, ρ_corrected_deph =
                    state_prep(ρ_mixed, angles_noisy_phase_real, angles_noisy_phase_imag, n_samples_populations, n_samples_phase_estimation_each)
                    # pure/dephased initial state sampling and phase correction

                mes_fidelity_pure[i] = fidelity(Ψ_mes, ρ_corrected_pure)
                mes_fidelity_deph[i] = fidelity(Ψ_mes, ρ_corrected_deph)
                    # true fidelities of the pure/dephased initial states after phase correction

                angles_noisy_single = noisy_angles_symmetric(angles_single, ϵ_angles)
                    # noisy angles for the single-setting DTQW

                pop_fs_coh_noisy_single_sampled_deph, fidelity_sampled_single_deph[i] = eval_single(
                    N, ρ_corrected_deph, pops_init_sampled_deph, angles_single, angles_noisy_single, n_samples_single_setting
                )
                pop_fs_coh_noisy_single_sampled_pure, fidelity_sampled_single_pure[i] = eval_single(
                    N, ρ_corrected_pure, pops_init_sampled_pure, angles_single, angles_noisy_single, n_samples_single_setting
                ) # single-setting DTQW evaluation and fidelity computation for pure/dephased initial states

                angles_noisy_compound = angles_compound(N, ϵ_angles)
                    # noisy angles for the compound-setting DTQWs

                pop_fs_coh_noisy_compound_sampled_deph = fs_pop_compound(
                    ρ_corrected_deph, angles_noisy_compound; n_samples = n_samples_compound_each
                )
                pop_fs_coh_noisy_compound_sampled_pure = fs_pop_compound(
                    ρ_corrected_pure, angles_noisy_compound; n_samples = n_samples_compound_each
                ) # calculate the sampled noisy final-state populations of the compound-setting DTQWs for pure/dephased initial states

                fidelity_sampled_compound_deph[i] = coherence_extraction_compound(
                    pops_init_sampled_deph, pop_fs_coh_noisy_compound_sampled_deph
                )
                fidelity_sampled_compound_pure[i] = coherence_extraction_compound(
                    pops_init_sampled_pure, pop_fs_coh_noisy_compound_sampled_pure
                ) # compound-setting DTQW evaluation and fidelity computation for pure/dephased initial states

                phase_measurement_pure[:, 1, i] .= Iterators.flatten(pops_fs_real_pure)
                phase_measurement_pure[:, 2, i] .= Iterators.flatten(pops_fs_imag_pure)
                phase_measurement_pure[:, 3, i] .= Iterators.flatten(pops_fs_real_sampled_pure)
                phase_measurement_pure[:, 4, i] .= Iterators.flatten(pops_fs_imag_sampled_pure)

                phase_measurement_deph[:, 1, i] .= Iterators.flatten(pops_fs_real_deph)
                phase_measurement_deph[:, 2, i] .= Iterators.flatten(pops_fs_imag_deph)
                phase_measurement_deph[:, 3, i] .= Iterators.flatten(pops_fs_real_sampled_deph)
                phase_measurement_deph[:, 4, i] .= Iterators.flatten(pops_fs_imag_sampled_deph)

                pop_init_sampled_pure_arr[:, i] = pops_init_sampled_pure
                pop_init_sampled_deph_arr[:, i] = pops_init_sampled_deph

                pop_fs_coh_noisy_single_sampled_deph_arr[i] = pop_fs_coh_noisy_single_sampled_deph
                pop_fs_coh_noisy_single_sampled_pure_arr[i] = pop_fs_coh_noisy_single_sampled_pure

                pop_fs_coh_noisy_compound_sampled_deph_arr[:, i] .= Iterators.flatten(pop_fs_coh_noisy_compound_sampled_deph)
                pop_fs_coh_noisy_compound_sampled_pure_arr[:, i] .= Iterators.flatten(pop_fs_coh_noisy_compound_sampled_pure)

                push!(angles_noisy_phase_real_arr, angles_noisy_phase_real)
                push!(angles_noisy_phase_imag_arr, angles_noisy_phase_imag)
                push!(angles_noisy_compound_arr, angles_noisy_compound)
                push!(angles_noisy_single_arr, angles_noisy_single)
            end
            output = mes_fidelity_deph, fidelity_sampled_single_deph, fidelity_sampled_compound_deph,
                mes_fidelity_pure, fidelity_sampled_single_pure, fidelity_sampled_compound_pure,
                phase_measurement_pure, phase_measurement_deph,
                pop_init_sampled_pure_arr, pop_init_sampled_deph_arr,
                pop_fs_coh_noisy_single_sampled_deph_arr, pop_fs_coh_noisy_compound_sampled_deph_arr,
                pop_fs_coh_noisy_single_sampled_pure_arr, pop_fs_coh_noisy_compound_sampled_pure_arr,
                angles_noisy_phase_real_arr, angles_noisy_phase_imag_arr,
                angles_noisy_compound_arr, angles_noisy_single_arr
            return output
        end
    end
   @time results = fetch.(tasks)

   # collect all the results

    mes_fidelity_deph = cat([results[i][1] for i in eachindex(data_chunks)]..., dims=1)
    fidelity_sampled_single_deph = cat([results[i][2] for i in eachindex(data_chunks)]..., dims=1)
    fidelity_sampled_compound_deph = cat([results[i][3] for i in eachindex(data_chunks)]..., dims=1)
    mes_fidelity_pure = cat([results[i][4] for i in eachindex(data_chunks)]..., dims=1)
    fidelity_sampled_single_pure = cat([results[i][5] for i in eachindex(data_chunks)]..., dims=1)
    fidelity_sampled_compound_pure = cat([results[i][6] for i in eachindex(data_chunks)]..., dims=1)
    phase_measurement_pure = cat([results[i][7] for i in eachindex(data_chunks)]..., dims=3)
    phase_measurement_deph = cat([results[i][8] for i in eachindex(data_chunks)]..., dims=3)
    pop_init_sampled_pure_arr = cat([results[i][9] for i in eachindex(data_chunks)]..., dims=2)
    pop_init_sampled_deph_arr = cat([results[i][10] for i in eachindex(data_chunks)]..., dims=2)
    pop_fs_coh_noisy_single_sampled_deph_arr = cat([results[i][11] for i in eachindex(data_chunks)]..., dims=1)
    pop_fs_coh_noisy_compound_sampled_deph_arr = cat([results[i][12] for i in eachindex(data_chunks)]..., dims=2)
    pop_fs_coh_noisy_single_sampled_pure_arr = cat([results[i][13] for i in eachindex(data_chunks)]..., dims=1)
    pop_fs_coh_noisy_compound_sampled_pure_arr = cat([results[i][14] for i in eachindex(data_chunks)]..., dims=2)
    angles_noisy_phase_real_arr = cat([results[i][15] for i in eachindex(data_chunks)]..., dims=1)
    angles_noisy_phase_imag_arr = cat([results[i][16] for i in eachindex(data_chunks)]..., dims=1)
    angles_noisy_compound_arr = cat([results[i][17] for i in eachindex(data_chunks)]..., dims=1)
    angles_noisy_single_arr = cat([results[i][18] for i in eachindex(data_chunks)]..., dims=1)

    data_path = string(@__DIR__,"/Data/N$(N)N_phase$(n_sample_phase).jld2")
    save(data_path,
        "phase_measurement_pure", phase_measurement_pure,
        "phase_measurement_deph", phase_measurement_deph,
        "angles_noisy_phase_real_arr", angles_noisy_phase_real_arr,
        "angles_noisy_phase_imag_arr", angles_noisy_phase_imag_arr,
        "angles_noisy_compound_arr", angles_noisy_compound_arr,
        "angles_noisy_single_arr", angles_noisy_single_arr,
        "pop_init_sampled_pure_arr", pop_init_sampled_pure_arr,
        "pop_init_sampled_deph_arr", pop_init_sampled_deph_arr,
        "pop_fs_coh_noisy_single_sampled_deph_arr", pop_fs_coh_noisy_single_sampled_deph_arr,
        "pop_fs_coh_noisy_single_sampled_pure_arr", pop_fs_coh_noisy_single_sampled_pure_arr,
        "pop_fs_coh_noisy_compound_sampled_deph_arr", pop_fs_coh_noisy_compound_sampled_deph_arr,
        "pop_fs_coh_noisy_compound_sampled_pure_arr", pop_fs_coh_noisy_compound_sampled_pure_arr,
        "mes_fidelity_deph", mes_fidelity_deph,
        "mes_fidelity_pure", mes_fidelity_pure,
        "fidelity_sampled_single_deph", fidelity_sampled_single_deph,
        "fidelity_sampled_single_pure", fidelity_sampled_single_pure,
        "fidelity_sampled_compound_deph", fidelity_sampled_compound_deph,
        "fidelity_sampled_compound_pure", fidelity_sampled_compound_pure,
    )

    F_med_single = median(fidelity_sampled_single_deph)
    F_std_single = std(fidelity_sampled_single_deph)
    F_med_compound = median(fidelity_sampled_compound_deph)
    F_std_compound = std(fidelity_sampled_compound_deph)
        # compute the median and standard deviation of the fidelities of single/compound-setting schemes for the dephased initial state

    F_med_single_pure = median(fidelity_sampled_single_pure)
    F_std_single_pure = std(fidelity_sampled_single_pure)
    F_med_compound_pure = median(fidelity_sampled_compound_pure)
    F_std_compound_pure = std(fidelity_sampled_compound_pure)
        # compute the median and standard deviation of the fidelities of single/compound-setting schemes for the pure initial state

    F_med_true = median(mes_fidelity_deph)
    F_std_true = std(mes_fidelity_deph)
        # compute the median and standard deviation of the true fidelities of the dephased initial state after phase correction

    F_med_true_pure = median(mes_fidelity_pure)
    F_std_true_pure = std(mes_fidelity_pure)
        # compute the median and standard deviation of the true fidelities of the pure initial state after phase correction

    println("")
    println("")
    println("N = $N | samples used: $(n_samples_populations + 2 * n_samples_phase_estimation_each + n_samples_single_setting)")
    println("#Samples Phase Estimation: $n_samples_phase_estimation_each")
    println("#Samples Populations : $n_samples_populations")
    println("#Samples Single Setting: $n_samples_single_setting")
    println("#Samples Compound Setting: $n_samples_compound_each")
    println("")
    println("Purity tr(ρ^2) = 1")
    println("True Fidelity ρ_corr to Ψ_MES = $(mes_fidelity_pure[1]) ± $F_std_true_pure | F_med = $F_med_true_pure")
    println("Single Setup: F = $(fidelity_sampled_single_pure[1]) ± $F_std_single_pure | F_med = $F_med_single_pure")
    println("Compound Setup: F = $(fidelity_sampled_compound_pure[1]) ± $F_std_compound_pure | F_med = $F_med_compound_pure")
    println("")
    println("Purity tr(ρ^2) ≈ 0.9")
    println("Max. Fidelity perfect phase correction = $(fidelity(Ψ_mes, ρ_perfect_phase_correction))")
    println("True Fidelity ρ_corr to Ψ_MES = $(mes_fidelity_deph[1]) ± $F_std_true | F_med = $F_med_true")
    println("Single Setup: F = $(fidelity_sampled_single_deph[1]) ± $F_std_single | F_med = $F_med_single")
    println("Compound Setup: F = $(fidelity_sampled_compound_deph[1]) ± $F_std_compound | F_med = $F_med_compound")

    flush(stdout)
end
