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

const n_samples_compound_each_arr = parse.(Int, split(chop(ARGS[7], head=1),','))
    # Number of samples for each of the (N-1) DTQWs in the compound setting

Random.seed!(8675309) # / Jenny dont change your number

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
ρ_pure = density_matrix(Ψ_mes)
angles_single = angles_single_setup(N)
    # angles for the single-setting DTQW
ρ_perfect_phase_correction = density_matrix_dephased(Ψ_mes, ϵ)
    # dephased density matrix, but with no relative phases, acting as a benchmark for the phase correction

for (idx_sample, n_samples_compound_each) in enumerate(n_samples_compound_each_arr)
    # iterate over the number of samples for each of the (N-1) DTQWs in the compound setting

    n_samples_populations = (N_samples - (N - 1) * n_samples_compound_each)
        # number of samples used for the initial-state populations

    n_samples_single_setting = (N - 1) * n_samples_compound_each
        # number of samples used for the one single-setting DTQW (total numer identical to the compound setting)

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

            local fidelity_sampled_single_pure = zeros(chunk_size)
            local fidelity_sampled_compound_pure = zeros(chunk_size)

            local pop_init_sampled_pure_arr = zeros(4 * N^2, chunk_size)
            local pop_init_sampled_deph_arr = zero(pop_init_sampled_pure_arr)

            local pop_fs_coh_noisy_single_sampled_deph_arr = zeros(chunk_size)
            local pop_fs_coh_noisy_single_sampled_pure_arr = zero(pop_fs_coh_noisy_single_sampled_deph_arr)

            local pop_fs_coh_noisy_compound_sampled_deph_arr = zeros(N - 1, chunk_size)
            local pop_fs_coh_noisy_compound_sampled_pure_arr = zero(pop_fs_coh_noisy_compound_sampled_deph_arr)

            local angles_noisy_compound_arr = []
            local angles_noisy_single_arr = []

            for i in 1:chunk_size
                # iterate over each independent experimental run
                pops_init_sampled_deph = populations(ρ_perfect_phase_correction, n_samples_populations)
                pops_init_sampled_pure = populations(ρ_pure, n_samples_populations)
                    # sample dephased and pure initial-state time-bin populations

                angles_noisy_single = noisy_angles_symmetric(angles_single, ϵ_angles)
                    # generate noisy angles for the single-setting DTQW

                pop_fs_coh_noisy_single_sampled_deph, fidelity_sampled_single_deph[i] = eval_single(
                    N, ρ_perfect_phase_correction, pops_init_sampled_deph, angles_single, angles_noisy_single, n_samples_single_setting
                )
                pop_fs_coh_noisy_single_sampled_pure, fidelity_sampled_single_pure[i] = eval_single(
                    N, ρ_pure, pops_init_sampled_pure, angles_single, angles_noisy_single, n_samples_single_setting
                )
                    # sampled dephased and pure final-state populations of the single-setting DTQW


                angles_noisy_compound = angles_compound(N, ϵ_angles)
                    # generate noisy angles for the compound-setting DTQW

                pop_fs_coh_noisy_compound_sampled_deph = fs_pop_compound(
                    ρ_perfect_phase_correction, angles_noisy_compound; n_samples = n_samples_compound_each
                )
                pop_fs_coh_noisy_compound_sampled_pure = fs_pop_compound(
                    ρ_pure, angles_noisy_compound; n_samples = n_samples_compound_each
                )
                    # calculate the dephased and pure noisy sampled final-state populations of the compound-setting DTQW

                fidelity_sampled_compound_deph[i] = coherence_extraction_compound(
                    pops_init_sampled_deph, pop_fs_coh_noisy_compound_sampled_deph
                )
                fidelity_sampled_compound_pure[i] = coherence_extraction_compound(
                    pops_init_sampled_pure, pop_fs_coh_noisy_compound_sampled_pure
                )
                    # compute the fidelity of the compound-setting DTQW for the dephased and pure initial states

                pop_init_sampled_pure_arr[:, i] = pops_init_sampled_pure
                pop_init_sampled_deph_arr[:, i] = pops_init_sampled_deph

                pop_fs_coh_noisy_single_sampled_deph_arr[i] = pop_fs_coh_noisy_single_sampled_deph
                pop_fs_coh_noisy_single_sampled_pure_arr[i] = pop_fs_coh_noisy_single_sampled_pure

                pop_fs_coh_noisy_compound_sampled_deph_arr[:, i] .= Iterators.flatten(pop_fs_coh_noisy_compound_sampled_deph)
                pop_fs_coh_noisy_compound_sampled_pure_arr[:, i] .= Iterators.flatten(pop_fs_coh_noisy_compound_sampled_pure)

                push!(angles_noisy_compound_arr, angles_noisy_compound)
                push!(angles_noisy_single_arr, angles_noisy_single)
            end
            output = fidelity_sampled_single_deph, fidelity_sampled_compound_deph,
                fidelity_sampled_single_pure, fidelity_sampled_compound_pure,
                pop_init_sampled_pure_arr, pop_init_sampled_deph_arr,
                pop_fs_coh_noisy_single_sampled_deph_arr, pop_fs_coh_noisy_compound_sampled_deph_arr,
                pop_fs_coh_noisy_single_sampled_pure_arr, pop_fs_coh_noisy_compound_sampled_pure_arr,
                angles_noisy_compound_arr, angles_noisy_single_arr
            return output
        end
    end
   @time results = fetch.(tasks)

    fidelity_sampled_single_deph = cat([results[i][1] for i in eachindex(data_chunks)]..., dims=1)
    fidelity_sampled_compound_deph = cat([results[i][2] for i in eachindex(data_chunks)]..., dims=1)
    fidelity_sampled_single_pure = cat([results[i][3] for i in eachindex(data_chunks)]..., dims=1)
    fidelity_sampled_compound_pure = cat([results[i][4] for i in eachindex(data_chunks)]..., dims=1)
    pop_init_sampled_pure_arr = cat([results[i][5] for i in eachindex(data_chunks)]..., dims=2)
    pop_init_sampled_deph_arr = cat([results[i][6] for i in eachindex(data_chunks)]..., dims=2)
    pop_fs_coh_noisy_single_sampled_deph_arr = cat([results[i][7] for i in eachindex(data_chunks)]..., dims=1)
    pop_fs_coh_noisy_compound_sampled_deph_arr = cat([results[i][8] for i in eachindex(data_chunks)]..., dims=2)
    pop_fs_coh_noisy_single_sampled_pure_arr = cat([results[i][9] for i in eachindex(data_chunks)]..., dims=1)
    pop_fs_coh_noisy_compound_sampled_pure_arr = cat([results[i][10] for i in eachindex(data_chunks)]..., dims=2)
    angles_noisy_compound_arr = cat([results[i][11] for i in eachindex(data_chunks)]..., dims=1)
    angles_noisy_single_arr = cat([results[i][12] for i in eachindex(data_chunks)]..., dims=1)

    #collect all results

    data_path = string(@__DIR__,"/Data/N$(N)N_comp$(n_samples_compound_each).jld2")
    save(data_path,
        "angles_noisy_compound_arr", angles_noisy_compound_arr,
        "angles_noisy_single_arr", angles_noisy_single_arr,
        "pop_init_sampled_pure_arr", pop_init_sampled_pure_arr,
        "pop_init_sampled_deph_arr", pop_init_sampled_deph_arr,
        "pop_fs_coh_noisy_single_sampled_deph_arr", pop_fs_coh_noisy_single_sampled_deph_arr,
        "pop_fs_coh_noisy_single_sampled_pure_arr", pop_fs_coh_noisy_single_sampled_pure_arr,
        "pop_fs_coh_noisy_compound_sampled_deph_arr", pop_fs_coh_noisy_compound_sampled_deph_arr,
        "pop_fs_coh_noisy_compound_sampled_pure_arr", pop_fs_coh_noisy_compound_sampled_pure_arr,
        "fidelity_sampled_single_deph", fidelity_sampled_single_deph,
        "fidelity_sampled_single_pure", fidelity_sampled_single_pure,
        "fidelity_sampled_compound_deph", fidelity_sampled_compound_deph,
        "fidelity_sampled_compound_pure", fidelity_sampled_compound_pure,
    )

    F_med_single = median(fidelity_sampled_single_deph)
    F_std_single = std(fidelity_sampled_single_deph)
    F_med_compound = median(fidelity_sampled_compound_deph)
    F_std_compound = std(fidelity_sampled_compound_deph)
        # compute the median and standard deviation of the fidelities of the single-setting and compound-setting schemes for the dephased initial state

    F_med_single_pure = median(fidelity_sampled_single_pure)
    F_std_single_pure = std(fidelity_sampled_single_pure)
    F_med_compound_pure = median(fidelity_sampled_compound_pure)
    F_std_compound_pure = std(fidelity_sampled_compound_pure)
        # compute the median and standard deviation of the fidelities of the single-setting and compound-setting schemes for the pure initial state


    println("")
    println("")
    println("N = $N | samples used: $(n_samples_populations + n_samples_single_setting)")
    println("#Samples Populations : $n_samples_populations")
    println("#Samples Single Setting: $n_samples_single_setting")
    println("#Samples Compound Setting: $n_samples_compound_each")
    println("")
    println("Purity tr(ρ^2) = 1")
    println("True Fidelity ρ_pure to Ψ_MES = $(fidelity(Ψ_mes, ρ_pure))")
    println("Single Setup: F = $(fidelity_sampled_single_pure[1]) ± $F_std_single_pure | F_med = $F_med_single_pure")
    println("Compound Setup: F = $(fidelity_sampled_compound_pure[1]) ± $F_std_compound_pure | F_med = $F_med_compound_pure")
    println("")
    println("Purity tr(ρ^2) ≈ 0.9")
    println("True Fidelity ρ_deph to Ψ_MES = $(fidelity(Ψ_mes, ρ_perfect_phase_correction))")
    println("Single Setup: F = $(fidelity_sampled_single_deph[1]) ± $F_std_single | F_med = $F_med_single")
    println("Compound Setup: F = $(fidelity_sampled_compound_deph[1]) ± $F_std_compound | F_med = $F_med_compound")

    flush(stdout)
end
