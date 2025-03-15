### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 537e52ec-cdb2-11ef-02ca-6fb2cd697476
begin
	using Revise
    import Pkg
    Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using LinearAlgebra, SparseArrays, TimeBinEncoding
	BLAS.set_num_threads(1)
end

# ╔═╡ f0d4a691-ca2d-4277-bf5c-ed19dda9a1ff
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2200px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
""" # adapt cell width and boundaries

# ╔═╡ 5c8d1f9d-520e-4149-8b66-1104379a1d95
N = 4 # Number of time bins

# ╔═╡ 63556294-c6f0-41f4-8ac0-76bcfebad954
begin
	d_local_hs_b = Int(N * (N + 1) / 2) # hilbert-space dimension of two identical photons in N time bins. No notion of short/long
	d_local_hs_bl = N * (2 * N + 1) # hilbert-space dimension of two identical photons in N time bins in full short/long notation
	d_full_hs_b = d_local_hs_b ^ 2 # hilbert-space dimension of 2 species of 2 photons each in N time bins. No notion of short/long
	d_full_hs_bl = d_local_hs_bl ^ 2 # hilbert-space dimension of 2 species of 2 photon each in N time binss in full short/long notion
end

# ╔═╡ 9f47bb00-ea2f-4c9a-a779-25a6ae5f3392
begin
	Ψ_MES = spzeros(ComplexF64, d_full_hs_bl) # maximally entangled state (MES)
	for l in 0:N - 1
		for m in l:N - 1
			j_super = lcmk2j_super_identical(N, l, 0, m, 0, l, 0, m, 0)
			Ψ_MES[j_super] = 1
				# all combinations of two time bins are equally probable
		end
	end
	normalize!(Ψ_MES) # normalization

	ρ_pure = density_matrix(Ψ_MES) # pure density matrix of the MES
	pops_pure = populations(ρ_pure) # populations of the MES
end

# ╔═╡ ade07eea-5d51-4517-902b-28be03648f84
begin
	ϵ = 0.05712 # dephasng strength
	ρ_mixed = density_matrix_dephased_krauss_identical(N, ρ_pure, ϵ)
		# dephased density matrix of two spdc pairs
	pops_mixed = populations(ρ_mixed)
		# populations in the tme-bin basis
end

# ╔═╡ f160d7e5-38b4-443e-9710-7dcf115e22d8
purity(ρ_mixed) # purity of the dephased state

# ╔═╡ 84610a30-6182-4083-8c26-01bc173b3b34
begin
	contr_j_idxs = correlated_short_bins_idxs_identical(N) # indices of contributing states to the fidelity
	contr_j_tuples=correlated_short_bins_tuples_identical(N; extract_diagonal=false)

	coh_pop = sum(pops_mixed[contr_j_idxs]) / d_local_hs_b

	projector_weigths_auto = combined_projector_weights_auto(N, [1, 1, 1, 1, 1])
	# weights of the different measured event propabilities in the compound quantity
	combined_weights_hom, pop_fs_hom = combined_weights_pops_hom(N, ρ_mixed, projector_weigths_auto)
		# combined weights and final state populations of the Hong-Ou-Mandel/compound scheme
end

# ╔═╡ 092f3a47-a100-4868-b3fc-25b396c7beaa
begin
	combined_weights_4b_simple, pop_fs_4b_simple = combined_weights_pops_4bins_all(N, ρ_mixed, projector_weigths_auto, phases=false)
		# compute the weights (decomposition) and the population (numerical value) of the simple four-bin-interference measurement scheme
	combined_weights_simple = combined_weights_hom + combined_weights_4b_simple
		# compute the decomposition of the final-state probability in terms of initial-state time-bin basis
	pop_fs_combined_simple = pop_fs_hom + pop_fs_4b_simple
		# compute the combined numerical value of the final-state probability
	visual_meas_coh_map_combined_identical(N, combined_weights_simple, contr_j_idxs)
		# Visualization of all cohereces contained in combined_weights_simple and of all coherences that get extracted
end

# ╔═╡ 983db4af-0f4c-49b5-9e4f-ea66a9d4aeca
begin
	combined_weights_4b_phases, pop_fs_4b_phases = combined_weights_pops_4bins_all(N, ρ_mixed, projector_weigths_auto, phases=true)
		# compute the weights (decomposition) and the population (numerical value) of the advanced four-bin-interference measurement scheme	with initial-state phase imprints
	combined_weights_phases = combined_weights_hom + combined_weights_4b_phases
		# compute the decomposition of the final-state probability in terms of initial-state time-bin basis
	pop_fs_combined_phases = pop_fs_hom + pop_fs_4b_phases
		# compute the combined numerical value of the final-state probability
	visual_meas_coh_map_combined_identical(N, combined_weights_phases, contr_j_idxs)
		# Visualization of all cohereces contained in combined_weights_phases and of all coherences that get extracted
end

# ╔═╡ 07a9c8c4-87f3-403c-9dfc-84686bfde09d
begin
	coh_extract_simple = combined_measurement_coherence_extraction_identical(
	    N,
	    combined_weights_simple,
	    pops_mixed,
	    pop_fs_combined_simple,
	    contr_j_tuples
	)# extract coherences four-bin measurements with no initial phase imprint and HOM measurements
	coh_extract_simple /= d_local_hs_b # normalization
	F_simple = coh_extract_simple + coh_pop # add contributions from populations and coherences
	println("Fidelity for simple method / no initial phase application: \nF = ", F_simple)
end

# ╔═╡ 5a377015-b2b0-4946-8833-fb47d4a94455
begin
	coh_extract_phases = combined_measurement_coherence_extraction_identical(
	    N,
	    combined_weights_phases,
	    pops_mixed,
	    pop_fs_combined_phases,
	    contr_j_tuples
	) # extract coherences four-bin measurements with differet initial phase imprints and HOM measurements
	coh_extract_phases /= d_local_hs_b # normalization
	F_phases = coh_extract_phases + coh_pop # add contributions from populations and coherences
	println("Fidelity for advanced method / with initial phase application: \nF = ", F_phases)
end

# ╔═╡ Cell order:
# ╠═537e52ec-cdb2-11ef-02ca-6fb2cd697476
# ╠═f0d4a691-ca2d-4277-bf5c-ed19dda9a1ff
# ╠═5c8d1f9d-520e-4149-8b66-1104379a1d95
# ╠═63556294-c6f0-41f4-8ac0-76bcfebad954
# ╠═9f47bb00-ea2f-4c9a-a779-25a6ae5f3392
# ╠═ade07eea-5d51-4517-902b-28be03648f84
# ╠═f160d7e5-38b4-443e-9710-7dcf115e22d8
# ╠═84610a30-6182-4083-8c26-01bc173b3b34
# ╠═092f3a47-a100-4868-b3fc-25b396c7beaa
# ╠═983db4af-0f4c-49b5-9e4f-ea66a9d4aeca
# ╠═07a9c8c4-87f3-403c-9dfc-84686bfde09d
# ╠═5a377015-b2b0-4946-8833-fb47d4a94455
