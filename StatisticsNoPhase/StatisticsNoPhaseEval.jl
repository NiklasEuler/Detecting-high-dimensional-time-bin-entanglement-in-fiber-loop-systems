### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ c451a096-ff32-4a2e-b083-1ba985c66100
begin
    import Pkg
    Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using StatsBase, FileIO, PyPlot, PyCall
	using TimeBinEncoding
end

# ╔═╡ 0142274a-9789-4cef-90d1-fe0fc928f8ee
@pyimport matplotlib.ticker as ticker

# ╔═╡ ab2169fe-c4d9-4581-8bbe-df96de940c8c
begin
	save_fig = true
    path = string(@__DIR__, "/Figs/")
	if save_fig
    	mkpath(path)
	end
end # make the path for saving the figure

# ╔═╡ 05533d08-b5ef-491c-93bd-648c67c19ae7
begin
	rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
	rcParams["font.size"] = 10
	rcParams["font.family"] = "serif"
	rcParams["font.serif"] = "STIX"
	rcParams["mathtext.fontset"] = "stix"
	rcParams["text.usetex"] = true
	rc("text.latex", preamble=raw"\usepackage{mathtools, bm, braket}")
	# load latex preamble

	cmap = plt.get_cmap("tab10")
	dimgrey = [70,105,105] / 255
	red = [255,0,0] / 255
	incr = red-dimgrey
end

# ╔═╡ 24c10bc2-2a06-4cc6-8ec7-68286957b60c
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2200px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ 873641be-e0cc-4be1-a810-fe6b21fd3e3b
begin
	n_resamples = 500 # number of reruns of the experiment
	N = 8 # number of time bins
	n_comp_arr = collect(140:60:560)
		# number of samples for each DTQW run in the compound scheme
	n_comp = length(n_comp_arr) # number of different sample distributions investigated
	N_samples = N * 1000 / 2
		# only half of the 8000 samples for populations and DTQW
end

# ╔═╡ c59eef98-3575-4663-80d4-24e1518568ca
begin
	n_DTQW = (N - 1) .* n_comp_arr
		# total number of samples for the DTQW settings for all sample distributions
	n_pop_arr = (N_samples .- n_DTQW)
		# number of samples for the population measurements for each sample distrbution
	sample_ratios = Dict(
		"DTQW"=>tuple(n_DTQW ./ N_samples),
		"Pop."=>tuple(n_pop_arr ./ N_samples),
	) # dictionary with the sample distributions
end 

# ╔═╡ 09559242-119c-4777-a90d-6256f54f700a
begin
	lower_threshold = 0.25
		# percantange of fidelity results indicated by lower errorbar
	upper_threshold = 0.75
		# percantange of fidelity results indicated by upper errorbar

	fidelity_sampled_single_pure = zeros(n_resamples, n_comp)
	fidelity_sampled_single_deph = zero(fidelity_sampled_single_pure)
	fidelity_sampled_comp_pure = zero(fidelity_sampled_single_pure)
	fidelity_sampled_comp_deph = zero(fidelity_sampled_single_pure)

	med_single_pure = zeros(n_comp)
	med_single_deph = zero(med_single_pure)
	med_comp_pure = zero(med_single_pure)
	med_comp_deph = zero(med_single_pure)

	lower_bound_single_pure = zero(med_single_pure)
	lower_bound_single_deph = zero(med_single_pure)
	lower_bound_comp_pure = zero(med_single_pure)
	lower_bound_comp_deph = zero(med_single_pure)

	upper_bound_single_pure = zero(med_single_pure)
	upper_bound_single_deph = zero(med_single_pure)
	upper_bound_comp_pure = zero(med_single_pure)
	upper_bound_comp_deph = zero(med_single_pure)

	for (idx_comp, N_samples_comp) in enumerate(n_comp_arr)
		# iterate over all sample distributions
		path = string(@__DIR__, "/Data/N$(N)N_comp$N_samples_comp.jld2")
		
		results = load(path) # load data
		fidelity_sampled_single_pure[:,idx_comp] = results["fidelity_sampled_single_pure"]
			# array of fidelity bound results obtained from a pure state using the single DTQW setup
		fidelity_sampled_single_deph[:,idx_comp] = results["fidelity_sampled_single_deph"]
			# array of fidelity bound results obtained from a dephased state using the single DTQW setup

		lower_quant_single_pure = quantile(fidelity_sampled_single_pure[:,idx_comp], lower_threshold)
		lower_quant_single_deph = quantile(fidelity_sampled_single_deph[:,idx_comp], lower_threshold)
			# compute the fidelity at the boundary of first and second quartile for the pure and dephased cases for the single DTQW setup
		
		higher_quant_single_pure = quantile(fidelity_sampled_single_pure[:,idx_comp], upper_threshold)
		higher_quant_single_deph = quantile(fidelity_sampled_single_deph[:,idx_comp], upper_threshold)
			# compute the fidelity at the boundary of the third and fourth quartile for the pure and dephased cases for the single DTQW setup

		med_single_pure[idx_comp] = median(fidelity_sampled_single_pure[:,idx_comp])
		med_single_deph[idx_comp] = median(fidelity_sampled_single_deph[:,idx_comp])
			# compute the corresponding fidelity medians

		upper_bound_single_pure[idx_comp] = higher_quant_single_pure - med_single_pure[idx_comp]
		upper_bound_single_deph[idx_comp] = higher_quant_single_deph - med_single_deph[idx_comp]
			# compute the deviation between the fidelity at the boundary of the third and fourth quartile and the median for the single DTQW setup

		lower_bound_single_pure[idx_comp] = med_single_pure[idx_comp] - lower_quant_single_pure
		lower_bound_single_deph[idx_comp] = med_single_deph[idx_comp] - lower_quant_single_deph
			# compute the deviation between the median and the fidelity at the boundary of the first and second quartile for the single DTQW setup

		fidelity_sampled_comp_pure[:,idx_comp] = results["fidelity_sampled_compound_pure"]
			# array of fidelity bound results obtained from a pure state using the compound DTQW setup
		fidelity_sampled_comp_deph[:,idx_comp] = results["fidelity_sampled_compound_deph"]
			# array of fidelity bound results obtained from a dephased state using the compound DTQW setup

		lower_quant_comp_pure = quantile(fidelity_sampled_comp_pure[:,idx_comp], lower_threshold)
		lower_quant_comp_deph = quantile(fidelity_sampled_comp_deph[:,idx_comp], lower_threshold)
			# compute the fidelity at the boundary of first and second quartile for the pure and dephased cases for the compound DTQW setup

		higher_quant_comp_pure = quantile(fidelity_sampled_comp_pure[:,idx_comp], upper_threshold)
		higher_quant_comp_deph = quantile(fidelity_sampled_comp_deph[:,idx_comp], upper_threshold)
			# compute the fidelity at the boundary of the third and fourth quartile for the pure and dephased cases for the compound DTQW setup

		med_comp_pure[idx_comp] = median(fidelity_sampled_comp_pure[:,idx_comp])
		med_comp_deph[idx_comp] = median(fidelity_sampled_comp_deph[:,idx_comp])
			# compute the corresponding fidelity medians

		upper_bound_comp_pure[idx_comp] = higher_quant_comp_pure - med_comp_pure[idx_comp]
		upper_bound_comp_deph[idx_comp] = higher_quant_comp_deph - med_comp_deph[idx_comp]
			# compute the deviation between the fidelity at the boundary of the third and fourth quartile and the median for the compound DTQW setup

		lower_bound_comp_pure[idx_comp] = med_comp_pure[idx_comp] - lower_quant_comp_pure
		lower_bound_comp_deph[idx_comp] = med_comp_deph[idx_comp] - lower_quant_comp_deph
			# compute the deviation between the median and the fidelity at the boundary of the first and second quartile for the compound DTQW setup

	end
	asymm_y_err_single_pure = [lower_bound_single_pure, upper_bound_single_pure]
	asymm_y_err_single_deph = [lower_bound_single_deph, upper_bound_single_deph]
	asymm_y_err_comp_pure = [lower_bound_comp_pure, upper_bound_comp_pure]
	asymm_y_err_comp_deph = [lower_bound_comp_deph, upper_bound_comp_deph]
		# asymmetric fidelity error bars for the fidelity bounds of the single-DTWQ scheme and the compound-DTQW scheme
end

# ╔═╡ dece04b6-e640-4890-977a-8d7df2594d6f
begin
	ϵ = 0.05215
	Ψ_mes = insert_initial_state(correlated_timebin_state(fill(1 / sqrt(N), N)))
		# maximally entangled state for N time bin
	ρ_pure = density_matrix(Ψ_mes)
		# corresponding pure density matrix
	ρ_perfect_phase_correction = density_matrix_dephased(Ψ_mes, ϵ)
		# corresponding dephased density matrix with zero initial phases
	fidel_deph = fidelity(Ψ_mes, ρ_perfect_phase_correction)
		# true fidelity between the dephased state and the MES
end

# ╔═╡ e454a30a-9d73-45fb-ad43-3d14c07fd427
fig2 = let
	column_width = 3.40457
	fig_height = column_width * 1.2
	cap_size = 3
	marker_size = 1.75
	lw_err = 1
	ls_err = ""
	fig, axs = subplots(3, figsize = (column_width, fig_height), sharex=true, height_ratios=[3,2,2])

	x_plot = 250:500:6000
	letter_arr = [L"$\bm{(a)}$", L"$\bm{(b)}$", L"$\bm{(c)}$"]

	for idx_ax in 1:2
		axs[idx_ax].grid()
		axs[idx_ax].text(-0.21, 0.96, letter_arr[idx_ax], transform=axs[idx_ax].transAxes)
		axs[idx_ax].yaxis.set_major_locator(ticker.MultipleLocator(0.05))
		axs[idx_ax].yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
		axs[idx_ax].set_ylabel(L"$F$")
		for lv in 7:7
			if idx_ax == 2
				continue
			end
			alpha = 0.9
			color = tuple([incr[j] * (lv - 1)/(N - 1) + dimgrey[j] for j in eachindex(incr)]...)
			axs[idx_ax].plot(x_plot, lv ./ fill(N, length(x_plot)), color = color, ls= (0, (3, 1)), lw= 0.75, alpha=alpha)
			# plot Bk thresholds
			y_pos = lv/N - 0.0175
			x_pos = 3900
			axs[idx_ax].text(x_pos, y_pos, latexstring("\$B_{$lv}\$"), color = color, fontsize="small", alpha=alpha)
		end
	end

	sing, cap, coll = axs[1].errorbar(n_DTQW, med_single_deph, yerr = asymm_y_err_single_deph, color = cmap(0), marker="<", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, label = "SS")
		# plot the fidelity bounds for the dephased state for the single-DTQW setup
	
	comp = axs[1].errorbar(n_DTQW, med_comp_deph, yerr = asymm_y_err_comp_deph, color = cmap(1), marker=">", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, label = "CS")
		# plot the fidelity bounds for the dephased state for the compound-DTQW setup

	corr, = axs[1].plot(x_plot, fill(fidel_deph, length(x_plot)), color = cmap(2), lw=lw_err, label = "true fidelity")
		# plot the true fidelity 

	axs[2].fill_between(n_DTQW, -lower_bound_single_deph, upper_bound_single_deph, zorder = 4)
	axs[2].fill_between(n_DTQW, -lower_bound_comp_deph, upper_bound_comp_deph, zorder = 4)
		# plot the 25%-75% confidence region for the single-DTQW and compound-DTQW fidelity bounds

	width = 125  # the width of the bars
	multiplier = -0.5

	keys = ["DTQW", "Pop."]
	colors = ["navy", "crimson", cmap(2)]

	for (key_idx, key) in enumerate(keys)
		value = sample_ratios[key][1]
		println(value)
		println(key)
	    offset = width * multiplier
	    rects = axs[3].bar(n_DTQW .+ offset, value, width, label=key, zorder = 8, color=colors[key_idx])
		# plot the percentage of the total sample number used for the different population and DTQW measurements 
	    multiplier += 1
	end

	axs[2].yaxis.set_major_locator(ticker.MultipleLocator(0.01))
	axs[2].yaxis.set_minor_locator(ticker.MultipleLocator(0.002))
	axs[2].set_ylabel(L"$Q_1\!-\!Q_3$"*"\n"*"envelope")

	axs[1].set_ylim([0.84, 0.96])
	axs[2].set_ylim([-0.011, 0.011])
	axs[3].set_ylim([0.0, 1])

	axs[3].grid()
	axs[3].text(-0.21, 1.2, letter_arr[3], transform=axs[3].transAxes)
	axs[3].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	axs[3].yaxis.set_minor_locator(ticker.MultipleLocator(0.04))
	axs[3].set_ylabel(L"$\mathrm{Sample~ratios}$")


	axs[1].set_xlim([n_DTQW[1] - 200, n_DTQW[end] + 200])
	axs[1].set_xticks(n_DTQW)
	axs[3].set_xlabel(L"$\#\mathrm{Samples~DTQW~total}$")

	axs[2].ticklabel_format(axis="y", style="sci", scilimits=(-2,-2))
	axs[3].ticklabel_format(axis="x", style="sci", scilimits=(1,1))

	axs[1].legend(ncol = 3, columnspacing = 1, loc = "lower center",
	 		  handlelength = 1.25, labelspacing = 0.3, framealpha = 1, borderpad = 0.35,
	 		  borderaxespad = 0., handletextpad = 0.2, fontsize = 8, bbox_to_anchor = (0.5,0.73), )
	axs[3].legend(ncol = 3, columnspacing = 1, loc = "lower center",
	 		  handlelength = 1.25, labelspacing = 0.3, framealpha = 1, borderpad = 0.35,
	 		  borderaxespad = 0., handletextpad = 0.2, fontsize = 8, bbox_to_anchor = (0.35, 0.75), )
	plt.tight_layout()
	plt.subplots_adjust(hspace=0.25)
	name = "FidelityEstimationStatistics_NoPhase"
	PyPlot.savefig(string(path, name, ".svg"), transparent=true)
	PyPlot.savefig(string(path, name, ".pdf"))
	PyPlot.close(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═c451a096-ff32-4a2e-b083-1ba985c66100
# ╠═0142274a-9789-4cef-90d1-fe0fc928f8ee
# ╠═ab2169fe-c4d9-4581-8bbe-df96de940c8c
# ╠═05533d08-b5ef-491c-93bd-648c67c19ae7
# ╠═24c10bc2-2a06-4cc6-8ec7-68286957b60c
# ╠═873641be-e0cc-4be1-a810-fe6b21fd3e3b
# ╠═c59eef98-3575-4663-80d4-24e1518568ca
# ╠═09559242-119c-4777-a90d-6256f54f700a
# ╠═dece04b6-e640-4890-977a-8d7df2594d6f
# ╠═e454a30a-9d73-45fb-ad43-3d14c07fd427
