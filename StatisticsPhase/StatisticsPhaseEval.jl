### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ b439b21e-b0c7-11ef-1116-5d992abc84fd
begin
    import Pkg
    Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using StatsBase, FileIO, PyPlot, PyCall
	using TimeBinEncoding
end

# ╔═╡ 920f20e1-6d84-4741-b3f4-fb97c1604043
@pyimport matplotlib.ticker as ticker

# ╔═╡ 7d236898-1f64-4b0b-9c70-82a6764a045f
begin
	save_fig = true
    path = string(@__DIR__, "/Figs/")
	if save_fig
    	mkpath(path)
	end
end # make the path for saving the figure

# ╔═╡ 6bcde7ef-3ce5-464c-91f6-be2635b98dc6
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

# ╔═╡ 784d33ae-850d-4e20-946e-945917337fb9
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

# ╔═╡ 4fbebd88-90a4-4a7b-bf85-f0eef1d9bb75
begin
	n_resamples = 500 # number of reruns of the experiment
	N = 8 # number of time bins
	N_samples = N * 1000 # total number of samples
end

# ╔═╡ 326fd4cf-aa14-4f07-acc7-d7fb451c9b0f
begin
	n_phase_arr = collect(200:350:2650)
		# number of samples used for each of the two DTQW settings needed for the phase estimation for each sample distribution
	n_phase_arr2 = 2 * n_phase_arr
		# total number of samples used for phase estimation for each sample distribution
	n_pop_arr = fill(2000, length(n_phase_arr))
		# number of samples used for population measurements for each sample distribution
	n_DTQW = (N_samples .- n_pop_arr .- n_phase_arr2)
		# total number of samples used for the DTQW measurements
	n_phase_sett = length(n_phase_arr)
	# number of differerent sample distribution investiagated
	sample_ratios = Dict(
		"DTQW"=>tuple(n_DTQW ./ N_samples),
		"Pop."=>tuple(n_pop_arr ./ N_samples),
	    "Phase"=>tuple(n_phase_arr2 ./ N_samples),
	) # dictionary with the sample distributions
end

# ╔═╡ f5030ce1-30f3-41b5-8b23-18f8433eefc2
begin
	lower_threshold = 0.25
		# percantange of fidelity results indicated by lower errorbar
	upper_threshold = 0.75
		# percantange of fidelity results indicated by upper errorbar
	
	fidelity_sampled_single_pure = zeros(n_resamples, n_phase_sett)
	fidelity_sampled_single_deph = zero(fidelity_sampled_single_pure)
	fidelity_sampled_comp_pure = zero(fidelity_sampled_single_pure)
	fidelity_sampled_comp_deph = zero(fidelity_sampled_single_pure)
	fidelity_true_corr_pure = zero(fidelity_sampled_single_pure)
	fidelity_true_corr_deph = zero(fidelity_sampled_single_pure)
	
	med_single_pure = zeros(n_phase_sett)
	med_single_deph = zero(med_single_pure)
	med_comp_pure = zero(med_single_pure)
	med_comp_deph = zero(med_single_pure)
	med_true_corr_pure = zero(med_single_pure)
	med_true_corr_deph = zero(med_single_pure)

	lower_bound_single_pure = zero(med_single_pure)
	lower_bound_single_deph = zero(med_single_pure)
	lower_bound_comp_pure = zero(med_single_pure)
	lower_bound_comp_deph = zero(med_single_pure)
	lower_bound_true_corr_pure = zero(med_single_pure)
	lower_bound_true_corr_deph = zero(med_single_pure)

	upper_bound_single_pure = zero(med_single_pure)
	upper_bound_single_deph = zero(med_single_pure)
	upper_bound_comp_pure = zero(med_single_pure)
	upper_bound_comp_deph = zero(med_single_pure)
	upper_bound_true_corr_pure = zero(med_single_pure)
	upper_bound_true_corr_deph = zero(med_single_pure)
	
	for (idx_phase, N_samples_phase) in enumerate(n_phase_arr)
		# iterate over all sample distributions, including phase estimation samples
		path = string(@__DIR__, "/Data/N$(N)N_phase$(N_samples_phase).jld2")
		
		results = load(path) #load data
		fidelity_sampled_single_pure[:,idx_phase] = results["fidelity_sampled_single_pure"]
			# array of fidelity bound results obtained from a pure state using the single DTQW setup and the phase-correction protocol
		fidelity_sampled_single_deph[:,idx_phase] = results["fidelity_sampled_single_deph"]
			# array of fidelity bound results obtained from a dephased state using the single DTQW setup and the phase-correction protocol
		
		lower_quant_single_pure = quantile(fidelity_sampled_single_pure[:,idx_phase], lower_threshold)
		lower_quant_single_deph = quantile(fidelity_sampled_single_deph[:,idx_phase], lower_threshold)
			# compute the fidelity at the boundary of first and second quartile for the pure and dephased cases for the single DTQW setup
		
		higher_quant_single_pure = quantile(fidelity_sampled_single_pure[:,idx_phase], upper_threshold)
		higher_quant_single_deph = quantile(fidelity_sampled_single_deph[:,idx_phase], upper_threshold)
			# compute the fidelity at the boundary of the third and fourth quartile for the pure and dephased cases for the single DTQW setup
	
		med_single_pure[idx_phase] = median(fidelity_sampled_single_pure[:,idx_phase])
		med_single_deph[idx_phase] = median(fidelity_sampled_single_deph[:,idx_phase])
			# compute the corresponding fidelity medians
	
		upper_bound_single_pure[idx_phase] = higher_quant_single_pure - med_single_pure[idx_phase]
		upper_bound_single_deph[idx_phase] = higher_quant_single_deph - med_single_deph[idx_phase]
			# compute the deviation between the fidelity at the boundary of the third and fourth quartile and the median for the single DTQW setup
	
		lower_bound_single_pure[idx_phase] = med_single_pure[idx_phase] - lower_quant_single_pure
		lower_bound_single_deph[idx_phase] = med_single_deph[idx_phase] - lower_quant_single_deph
			# compute the deviation between the median and the fidelity at the boundary of the first and second quartile for the single DTQW setup
	
		fidelity_sampled_comp_pure[:,idx_phase] = results["fidelity_sampled_compound_pure"]
			# array of fidelity bound results obtained from a pure state using the compound DTQW setup and the phase-correction protocol
		fidelity_sampled_comp_deph[:,idx_phase] = results["fidelity_sampled_compound_deph"]
			# array of fidelity bound results obtained from a dephased state using the compound DTQW setup and the phase-correction protocol
	
		lower_quant_comp_pure = quantile(fidelity_sampled_comp_pure[:,idx_phase], lower_threshold)
		lower_quant_comp_deph = quantile(fidelity_sampled_comp_deph[:,idx_phase], lower_threshold)
			# compute the fidelity at the boundary of first and second quartile for the pure and dephased cases for the compound DTQW setup
		
		higher_quant_comp_pure = quantile(fidelity_sampled_comp_pure[:,idx_phase], upper_threshold)
		higher_quant_comp_deph = quantile(fidelity_sampled_comp_deph[:,idx_phase], upper_threshold)
			# compute the fidelity at the boundary of the third and fourth quartile for the pure and dephased cases for the compound DTQW setup
	
		med_comp_pure[idx_phase] = median(fidelity_sampled_comp_pure[:,idx_phase])
		med_comp_deph[idx_phase] = median(fidelity_sampled_comp_deph[:,idx_phase])
			# compute the corresponding fidelity medians
	
		upper_bound_comp_pure[idx_phase] = higher_quant_comp_pure - med_comp_pure[idx_phase]
		upper_bound_comp_deph[idx_phase] = higher_quant_comp_deph - med_comp_deph[idx_phase]
			# compute the deviation between the fidelity at the boundary of the third and fourth quartile and the median for the compound DTQW setup
	
		lower_bound_comp_pure[idx_phase] = med_comp_pure[idx_phase] - lower_quant_comp_pure
		lower_bound_comp_deph[idx_phase] = med_comp_deph[idx_phase] - lower_quant_comp_deph
			# compute the deviation between the median and the fidelity at the boundary of the first and second quartile for the compound DTQW setup

		fidelity_true_corr_pure[:,idx_phase] = results["mes_fidelity_pure"]
			# array of true fidelities obtained from a pure random-phase initial state with the phase-free MES after phase correction
		fidelity_true_corr_deph[:,idx_phase] = results["mes_fidelity_deph"]
			# array of true fidelities obtained from a dephased random-phase initial state with the phase-free MES after phase correction
	
		lower_quant_corr_pure = quantile(fidelity_true_corr_pure[:,idx_phase], lower_threshold)
		lower_quant_corr_deph = quantile(fidelity_true_corr_deph[:,idx_phase], lower_threshold)
			# compute the fidelity at the boundary of first and second quartile for the pure and dephased cases for the true phase-corrected fidelity
			
		higher_quant_corr_pure = quantile(fidelity_true_corr_pure[:,idx_phase], upper_threshold)
		higher_quant_corr_deph = quantile(fidelity_true_corr_deph[:,idx_phase], upper_threshold)
			# compute the fidelity at the boundary of the third and fourth quartile for the pure and dephased cases for the true phase-corrected fidelity
	
		med_true_corr_pure[idx_phase] = median(fidelity_true_corr_pure[:,idx_phase])
		med_true_corr_deph[idx_phase] = median(fidelity_true_corr_deph[:,idx_phase])
			# compute the corresponding fidelity medians
	
		upper_bound_true_corr_pure[idx_phase] = higher_quant_corr_pure - med_true_corr_pure[idx_phase]
		upper_bound_true_corr_deph[idx_phase] = higher_quant_corr_deph - med_true_corr_deph[idx_phase]
			# compute the deviation between the fidelity at the boundary of the third and fourth quartile and the median for the true phase-corrected fidelity
	
		lower_bound_true_corr_pure[idx_phase] = med_true_corr_pure[idx_phase] - lower_quant_corr_pure
		lower_bound_true_corr_deph[idx_phase] = med_true_corr_deph[idx_phase] - lower_quant_corr_deph
			# compute the deviation between the median and the fidelity at the boundary of the first and second quartile for the true phase-corrected fidelity
	end
	asymm_y_err_single_pure = [lower_bound_single_pure, upper_bound_single_pure]
	asymm_y_err_single_deph = [lower_bound_single_deph, upper_bound_single_deph]
	asymm_y_err_comp_pure = [lower_bound_comp_pure, upper_bound_comp_pure]
	asymm_y_err_comp_deph = [lower_bound_comp_deph, upper_bound_comp_deph]
	asymm_y_err_corr_pure = [lower_bound_true_corr_pure, upper_bound_true_corr_pure]
	asymm_y_err_corr_deph = [lower_bound_true_corr_deph, upper_bound_true_corr_deph]
		# asymmetric fidelity error bars for the fidelity bounds of the single-DTWQ scheme, the compound-DTQW scheme, an the true phase-corrected fidelities
end

# ╔═╡ 10afd567-694d-49d3-86ca-dfcb314d15b0
begin
	ϵ = 0.05215
	Ψ_mes = insert_initial_state(correlated_timebin_state(fill(1 / sqrt(N), N)))
	ρ_pure = density_matrix(Ψ_mes)
	ρ_perfect_phase_correction = density_matrix_dephased(Ψ_mes, ϵ)
	fidel_deph = fidelity(Ψ_mes, ρ_perfect_phase_correction)
end

# ╔═╡ df87750e-b7da-4fd8-b6d6-31a6de9e072b
fig2 = let
	column_width = 3.40457
	fig_height = column_width * 1.2
	cap_size = 3
	marker_size = 1.75
	lw_err = 1
	ls_err = ""
	fig, axs = subplots(3, figsize = (column_width, fig_height), sharex=true, height_ratios=[3,2,2])

	x_plot = -100:500:7000
	letter_arr = [L"$\bm{(a)}$", L"$\bm{(b)}$", L"$\bm{(c)}$"]
			
	for idx_ax in 1:2
		axs[idx_ax].grid()
		axs[idx_ax].text(-0.21, 0.96, letter_arr[idx_ax], transform=axs[idx_ax].transAxes)
		axs[idx_ax].yaxis.set_major_locator(ticker.MultipleLocator(0.1))
		axs[idx_ax].yaxis.set_minor_locator(ticker.MultipleLocator(0.02))
		axs[idx_ax].set_ylabel(L"$F$")
		for lv in 5:7
			if idx_ax == 2
				continue
			end
			alpha = 0.9
			color = tuple([incr[j] * (lv - 1)/(N - 1) + dimgrey[j] for j in eachindex(incr)]...)
			axs[idx_ax].plot(x_plot, lv ./ fill(N, length(x_plot)), color = color, ls= (0, (3, 1)), lw= 0.75, alpha=alpha)
			y_pos = lv/N - 0.04
			x_pos = 5500
			if lv < 5
				continue
			end
			axs[idx_ax].text(x_pos, y_pos, latexstring("\$B_{$lv}\$"), color = color, fontsize="small", alpha=alpha)
		end
	end
	corr, = axs[1].errorbar(n_phase_arr2, med_true_corr_deph, yerr = asymm_y_err_corr_deph, color = cmap(2), marker="v", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, label = "true fidelity")
		# plot the dephased true phase-corrected fidelity medians
	
	sing, = axs[1].errorbar(n_phase_arr2, med_single_deph, yerr = asymm_y_err_single_deph, color = cmap(0), marker="<", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, label = "SS")
		# plot the fidelity-bound medians for the dephased state for the single-DTQW setup including phase correction

	comp, = axs[1].errorbar(n_phase_arr2, med_comp_deph, yerr = asymm_y_err_comp_deph, color = cmap(1), marker=">", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, label = "CS")
		# plot the fidelity-bound medians for the dephased state for the compound-DTQW setup including phase correction
	
	corr, = axs[1].plot(x_plot, fill(fidel_deph, length(x_plot)), color = cmap(2), lw=lw_err, label = "true fidelity\n(phase free)")
		# plot the dephased true fidelities for all N with no initial phases

	axs[2].fill_between(n_phase_arr2, -lower_bound_true_corr_deph, upper_bound_true_corr_deph, alpha=0.8, facecolor=cmap(2), edgecolor=cmap(2), lw=0.75, zorder=4)
		# plot the 25%-75% confidence region for the single-DTQW and compound-DTQW fidelity bounds
	axs[2].set_ylabel(L"$Q_1\!-\!Q_3$"*"\n"*"envelope")

	
	width = 200  # the width of the bars
	multiplier = -1

	keys = ["DTQW", "Pop.", "Phase"]
	colors = ["navy", "crimson", cmap(2)]
	
	for (key_idx, key) in enumerate(keys)
		value = sample_ratios[key][1]
		println(value)
		println(key)
	    offset = width * multiplier
	    rects = axs[3].bar(n_phase_arr2 .+ offset, value, width, label=key, zorder = 8, color=colors[key_idx])
	    multiplier += 1
	end
	
	axs[1].set_ylim([0.54, 0.96])
	axs[2].set_ylim([-0.15, 0.11])
	axs[3].set_ylim([0.0, 0.8])

	axs[3].grid()
	axs[3].text(-0.21, 1.2, letter_arr[3], transform=axs[3].transAxes)
	axs[3].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	axs[3].yaxis.set_minor_locator(ticker.MultipleLocator(0.04))
	axs[3].set_ylabel(L"$\mathrm{Sample~ratios}$")

	
	axs[1].set_xlim([n_phase_arr2[1] - 500, n_phase_arr2[end] + 600])
	axs[1].set_xticks(n_phase_arr2)
	axs[3].set_xlabel(L"$\#\mathrm{Samples~phase~est.~total}$")#, labelpad = 0)


	axs[2].ticklabel_format(axis="y", style="sci", scilimits=(-1,-1))	
	axs[3].ticklabel_format(axis="x", style="sci", scilimits=(2,2))
	
	axs[1].legend(ncol = 3, columnspacing = 0.75, loc = "lower center",
	 		  handlelength = 1.25, labelspacing = 0.5, framealpha = 1, borderpad = 0.35,
	 		  borderaxespad = 0., handletextpad = 0.2, fontsize = 8, bbox_to_anchor = (0.5,0.025), )
	axs[3].legend(ncol = 3, columnspacing = 1, loc = "upper center",
	 		  handlelength = 1.25, labelspacing = 0.5, framealpha = 1, borderpad = 0.35,
	 		  borderaxespad = 0., handletextpad = 0.2, fontsize = 8, bbox_to_anchor = (0.54,0.975), )
	plt.tight_layout()
	plt.subplots_adjust(hspace=0.25)
	name = "FidelityEstimationStatistics_Phase"
	PyPlot.savefig(string(path, name, ".svg"), transparent=true)
	PyPlot.savefig(string(path, name, ".pdf"))
	PyPlot.close(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═b439b21e-b0c7-11ef-1116-5d992abc84fd
# ╠═920f20e1-6d84-4741-b3f4-fb97c1604043
# ╠═7d236898-1f64-4b0b-9c70-82a6764a045f
# ╠═6bcde7ef-3ce5-464c-91f6-be2635b98dc6
# ╠═784d33ae-850d-4e20-946e-945917337fb9
# ╠═4fbebd88-90a4-4a7b-bf85-f0eef1d9bb75
# ╠═326fd4cf-aa14-4f07-acc7-d7fb451c9b0f
# ╠═f5030ce1-30f3-41b5-8b23-18f8433eefc2
# ╠═10afd567-694d-49d3-86ca-dfcb314d15b0
# ╠═df87750e-b7da-4fd8-b6d6-31a6de9e072b
