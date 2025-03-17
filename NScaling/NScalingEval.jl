### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ b7bf6320-ae60-11ef-06f3-2d4cb0982e2d
begin
    import Pkg
    Pkg.activate(Base.current_project())
	Pkg.instantiate()
	using StatsBase, FileIO, PyPlot, PyCall
	using TimeBinEncoding
end

# ╔═╡ 38280dd2-5a30-46a8-8f08-a3aef31dc081
begin
	@pyimport matplotlib.ticker as ticker
	@pyimport matplotlib.patches as patches
end

# ╔═╡ a9fd6f6a-9461-4b4d-8e11-0c64d7e1918d
begin
	save_fig = true
    path = string(@__DIR__, "/Figs/")
	if save_fig
    	mkpath(path)
	end
end # make the path for saving the figure

# ╔═╡ a9bdea14-43a0-4550-bee3-4691141e19e3
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

# ╔═╡ 46ad260d-44ae-4eb3-afb9-2a16e72e2e2f
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2200px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
""" # adapt cell widths and padding

# ╔═╡ 1f8b29b6-03f0-46fa-a455-ca789a4725d0
begin
	n_resamples = 500 # number of reruns of the experiment
	N_arr = [2, 4, 8, 16] # number of time bins
	n_N = length(N_arr) # number of different values of N investigated
	N_max = N_arr[end] # largest N investigated
end

# ╔═╡ 537b1e3e-d4e7-4bd3-872d-83b9fbd3c15f
begin
	ϵ_arr = [0.06905, 0.054836, 0.05215, 0.05152]
		# array of epsilon values, gauged such that tr(ρ^2) ≈ 0.9 for all N
	fidel_deph_true = zero(ϵ_arr)
	for (idx_N, N) in enumerate(N_arr)
		Ψ_mes = insert_initial_state(correlated_timebin_state(fill(1 / sqrt(N), N)))
			# maximally entangled state for N time bin
		ρ_pure = density_matrix(Ψ_mes)
			# corresponding pure density matrix
		ρ_perfect_phase_correction = density_matrix_dephased(Ψ_mes, ϵ_arr[idx_N])
			# corresponding dephased density matrix with zero initial phases
		fidel_deph_true[idx_N] = fidelity(Ψ_mes, ρ_perfect_phase_correction)
			# fidelity between the pure and dephased states for all N
	end
end

# ╔═╡ 3d043162-1554-4ef1-9709-73bfacc31346
begin
	lower_threshold = 0.25
		# percantange of fidelity results indicated by lower errorbar
	upper_threshold = 0.75
		# percantange of fidelity results indicated by upper errorbar


	fidelity_sampled_single_pure = zeros(n_resamples, n_N)
	fidelity_sampled_single_deph = zero(fidelity_sampled_single_pure)
	fidelity_sampled_comp_pure = zero(fidelity_sampled_single_pure)
	fidelity_sampled_comp_deph = zero(fidelity_sampled_single_pure)
	fidelity_true_corr_pure = zero(fidelity_sampled_single_pure)
	fidelity_true_corr_deph = zero(fidelity_sampled_single_pure)

	med_single_pure = zeros(n_N)
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

	for (idx_N,N) in enumerate(N_arr)
		# iterate over all N ∈ [2,4,8,16]
		path = string(@__DIR__, "/Data/N$N.jld2")
		results = load(path) # load data
		fidelity_sampled_single_pure[:,idx_N] = results["fidelity_sampled_single_pure"]
			# array of fidelity bound results obtained from a pure state using the single DTQW setup
		fidelity_sampled_single_deph[:,idx_N] = results["fidelity_sampled_single_deph"]
			# array of fidelity bound results obtained from a dephased state using the single DTQW setup

		lower_quant_single_pure = quantile(fidelity_sampled_single_pure[:,idx_N], lower_threshold)
		lower_quant_single_deph = quantile(fidelity_sampled_single_deph[:,idx_N], lower_threshold)
			# compute the fidelity at the boundary of first and second quartile for the pure and dephased cases for the single DTQW setup

		higher_quant_single_pure = quantile(fidelity_sampled_single_pure[:,idx_N], upper_threshold)
		higher_quant_single_deph = quantile(fidelity_sampled_single_deph[:,idx_N], upper_threshold)
			# compute the fidelity at the boundary of the third and fourth quartile for the pure and dephased cases for the single DTQW setup

		med_single_pure[idx_N] = median(fidelity_sampled_single_pure[:,idx_N])
		med_single_deph[idx_N] = median(fidelity_sampled_single_deph[:,idx_N])
			# compute the corresponding fidelity medians

		upper_bound_single_pure[idx_N] = higher_quant_single_pure - med_single_pure[idx_N]
		upper_bound_single_deph[idx_N] = higher_quant_single_deph - med_single_deph[idx_N]
			# compute the deviation between the fidelity at the boundary of the third and fourth quartile and the median for the single DTQW setup

		lower_bound_single_pure[idx_N] = med_single_pure[idx_N] - lower_quant_single_pure
		lower_bound_single_deph[idx_N] = med_single_deph[idx_N] - lower_quant_single_deph
			# compute the deviation between the median and the fidelity at the boundary of the first and second quartile for the single DTQW setup

		fidelity_sampled_comp_pure[:,idx_N] = results["fidelity_sampled_compound_pure"]
			# array of fidelity bound results obtained from a pure state using the compound DTQW setup
		fidelity_sampled_comp_deph[:,idx_N] = results["fidelity_sampled_compound_deph"]
			# array of fidelity bound results obtained from a dephased state using the compound DTQW setup

		lower_quant_comp_pure = quantile(fidelity_sampled_comp_pure[:,idx_N], lower_threshold)
		lower_quant_comp_deph = quantile(fidelity_sampled_comp_deph[:,idx_N], lower_threshold)
			# compute the fidelity at the boundary of first and second quartile for the pure and dephased cases for the compound DTQW setup

		higher_quant_comp_pure = quantile(fidelity_sampled_comp_pure[:,idx_N], upper_threshold)
		higher_quant_comp_deph = quantile(fidelity_sampled_comp_deph[:,idx_N], upper_threshold)
			# compute the fidelity at the boundary of the third and fourth quartile for the pure and dephased cases for the compound DTQW setup

		med_comp_pure[idx_N] = median(fidelity_sampled_comp_pure[:,idx_N])
		med_comp_deph[idx_N] = median(fidelity_sampled_comp_deph[:,idx_N])
			# compute the corresponding fidelity medians

		upper_bound_comp_pure[idx_N] = higher_quant_comp_pure - med_comp_pure[idx_N]
		upper_bound_comp_deph[idx_N] = higher_quant_comp_deph - med_comp_deph[idx_N]
			# compute the deviation between the fidelity at the boundary of the third and fourth quartile and the median for the compound DTQW setup

		lower_bound_comp_pure[idx_N] = med_comp_pure[idx_N] - lower_quant_comp_pure
		lower_bound_comp_deph[idx_N] = med_comp_deph[idx_N] - lower_quant_comp_deph
			# compute the deviation between the median and the fidelity at the boundary of the first and second quartile for the compound DTQW setup

		fidelity_true_corr_pure[:,idx_N] = results["mes_fidelity_pure"]
			# array of true fidelities obtained from a pure random-phase initial state with the phase-free MES after phase correction
		fidelity_true_corr_deph[:,idx_N] = results["mes_fidelity_deph"]
			# array of true fidelities obtained from a dephased random-phase initial state with the phase-free MES after phase correction

		lower_quant_corr_pure = quantile(fidelity_true_corr_pure[:,idx_N], lower_threshold)
		lower_quant_corr_deph = quantile(fidelity_true_corr_deph[:,idx_N], lower_threshold)
			# compute the fidelity at the boundary of first and second quartile for the pure and dephased cases for the true phase-corrected fidelity

		higher_quant_corr_pure = quantile(fidelity_true_corr_pure[:,idx_N], upper_threshold)
		higher_quant_corr_deph = quantile(fidelity_true_corr_deph[:,idx_N], upper_threshold)
			# compute the fidelity at the boundary of the third and fourth quartile for the pure and dephased cases for the true phase-corrected fidelity

		med_true_corr_pure[idx_N] = median(fidelity_true_corr_pure[:,idx_N])
		med_true_corr_deph[idx_N] = median(fidelity_true_corr_deph[:,idx_N])
			# compute the corresponding fidelity medians

		upper_bound_true_corr_pure[idx_N] = higher_quant_corr_pure - med_true_corr_pure[idx_N]
		upper_bound_true_corr_deph[idx_N] = higher_quant_corr_deph - med_true_corr_deph[idx_N]
			# compute the deviation between the fidelity at the boundary of the third and fourth quartile and the median for the true phase-corrected fidelity

		lower_bound_true_corr_pure[idx_N] = med_true_corr_pure[idx_N] - lower_quant_corr_pure
		lower_bound_true_corr_deph[idx_N] = med_true_corr_deph[idx_N] - lower_quant_corr_deph
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

# ╔═╡ e1027453-7355-4ceb-8cbf-6fbcb89e629b
fig = begin
	column_width = 3.40457
	fig_height = column_width * 1.2
	cap_size = 3
	marker_size = 2
	lw_err = 1
	ls_err = ""
	fig, axs = subplots(2, figsize = (column_width, fig_height), sharex=true)

	x_plot = 1:0.1:4.5
	letter_arr = [L"$\bm{(a)}$", L"$\bm{(b)}$"]

	xN = 1:4

	for idx_ax in 1:2
		axs[idx_ax].set_ylim([0.7, 1.0175])
		axs[idx_ax].grid()
		axs[idx_ax].set_ylabel(L"$F$")
		axs[idx_ax].yaxis.set_major_locator(ticker.MultipleLocator(0.1))
		axs[idx_ax].yaxis.set_minor_locator(ticker.MultipleLocator(0.02))
		axs[idx_ax].text(0.1, 0.98, letter_arr[idx_ax])
		for lv in 3:N_max
			if lv % 3 == 0
				alpha = 0.9
			else
				alpha = 0.25
			end
			color = tuple([incr[j] * (lv - 1)/(N_max - 1) + dimgrey[j] for j in eachindex(incr)]...)
			axs[idx_ax].plot(x_plot, lv ./ (2 .^ (x_plot)), color = color, ls= (0, (3, 1)), lw= 0.75, alpha=alpha)
				# plot the Bk thresholds
			if lv % 3 == 0
				if 3 < lv < 12
					y_pos = 0.85
					x_pos = log2(lv/y_pos) -  0.21
				elseif lv < 10
					y_pos = 0.975
					x_pos = log2(lv/y_pos) -  0.215
				else
					y_pos = 0.975
					x_pos = log2(lv/y_pos) - 0.28
				end
				axs[idx_ax].text(x_pos, y_pos, latexstring("\$B_{$lv}\$"), color = color, fontsize="small", alpha=alpha)
			end
		end
	end
	axs[1].scatter(xN, fill(1, length(xN)), marker = "x", lw=lw_err, label = "true fidelity\n (phase free)", zorder=2, color="darkgreen")
	axs[2].scatter(xN, fidel_deph_true, marker = "x", lw=lw_err, label = "true fidelity\n (phase free)", zorder=2, color="darkgreen")
		# plot the pure and dephased true fidelities of a phase-free initial state

	axs[1].errorbar(xN, med_true_corr_pure, yerr = asymm_y_err_corr_pure, color = cmap(2), marker="v", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, label = "true fidelity")
	axs[2].errorbar(xN, med_true_corr_deph, yerr = asymm_y_err_corr_deph, color = cmap(2), marker="v", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, label = "true fidelity")
		# plot the pure and dephased true fidelity medians after application of the phase-correction protocol
	axs[1].errorbar(xN, med_single_pure, yerr = asymm_y_err_single_pure, color = cmap(0), marker="<", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, label = "SS")
	axs[2].errorbar(xN, med_single_deph, yerr = asymm_y_err_single_deph, color = cmap(0), marker="<", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, label = "SS")
		# plot the true and dephased fidelity bound-medians from the single-DTQW setup, including phase correction
	axs[1].errorbar(xN, med_comp_pure, yerr = asymm_y_err_comp_pure, color = cmap(1), marker=">", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, label = "CS")
	axs[2].errorbar(xN, med_comp_deph, yerr = asymm_y_err_comp_deph, color = cmap(1), marker=">", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, label = "CS")
		# plot the true and dephased fidelity bound-medians from the compound-DTQW setup, including phase correction

	axs[1].set_xlim([0.8, 4.15])
	axs[1].set_xticks(1:n_N, [2^n for n in 1:n_N])
	axs[2].set_xlabel(L"$N$", labelpad = 0)

	x1, x2, y1, y2 = 0.9, 1.1, 0.97, 1.0125
	ix1, iy1, ixw, iyh = 0.9, 0.715, 0.5, 0.17
	axins1 = axs[1].inset_axes(
    	[ix1, iy1, ixw, iyh], xlim=(x1, x2), xticklabels=[], transform=axs[1].transData
	)
	axins1.set_ylim([0.97, 1.003])
	axins1.yaxis.set_major_locator(ticker.MultipleLocator(0.01))
	axins1.yaxis.set_minor_locator(ticker.MultipleLocator(0.005))
	axins1.yaxis.tick_right()
	yticks = [0.97, 0.98, 0.99, 1.0]
	axins1.set_yticks(yticks)
	axins1.set_yticklabels(yticks, fontsize=8)
	axins1.tick_params(axis="y", which="major", pad=2)
	axins1.grid()
	axins1.set_xticks([], [])

	zoom_axs1 = patches.Rectangle((x1, y1), x2-x1, y2-y1, fill=false, lw=lw_err, color = "grey")
	axs[1].plot([x1, ix1],[y1, iy1 + iyh], color="grey", lw=1)
	axs[1].plot([x2, ix1 + ixw],[y2, iy1 + iyh], color="grey", lw=1)


	idx_N = 1
	axins1.errorbar(xN, med_single_pure, yerr = asymm_y_err_single_pure, color = cmap(0), marker="<", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, zorder=3)
	axins1.errorbar(xN, med_comp_pure, yerr = asymm_y_err_comp_pure, color = cmap(1), marker=">", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, zorder=3)
	axins1.errorbar(xN, med_true_corr_pure, yerr = asymm_y_err_corr_pure, color = cmap(2), marker="v", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, zorder=3)
	axins1.scatter(xN, fill(1, length(xN)), marker = "x", lw=lw_err, zorder=2, color="darkgreen")
		# zoom in on N=2 data for a pure initial state

	x1, x2, y1, y2 = 0.9, 1.1, 0.9, 0.96  # subregion of the first ax
	ix1, iy1, ixw, iyh =  0.9, 0.71, 0.5, 0.17
	axins2 = axs[2].inset_axes(
    	[ix1, iy1, ixw, iyh], xlim=(x1, x2), xticklabels=[], transform=axs[2].transData
	)
	axins2.set_ylim([0.8975, 0.9525])
	axins2.yaxis.set_major_locator(ticker.MultipleLocator(0.01))
	axins2.yaxis.set_minor_locator(ticker.MultipleLocator(0.005))
	axins2.yaxis.tick_right()
	yticks = 0.9:0.02:0.95
	axins2.set_yticks(yticks)
	axins2.set_yticklabels(yticks, fontsize=8)
	axins2.set_xticks([], [])

	axins2.tick_params(axis="x", which="both", bottom="false", top="false", labelbottom="false")

	axins2.grid()

	zoom_axs2 = patches.Rectangle((x1, y1), x2-x1, y2-y1, fill=false, lw=lw_err, color = "grey")
	axs[2].plot([x1, ix1],[y1, iy1 + iyh], color="grey", lw=1)
	axs[2].plot([x2, ix1 + ixw],[y2, iy1 + iyh], color="grey", lw=1)

	axins2.errorbar(xN, med_single_deph, yerr = asymm_y_err_single_deph, color = cmap(0), marker="<", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, zorder=3)
	axins2.errorbar(xN, med_comp_deph, yerr = asymm_y_err_comp_deph, color = cmap(1), marker=">", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, zorder=3)
	axins2.errorbar(xN, med_true_corr_deph, yerr = asymm_y_err_corr_deph, color = cmap(2), marker="v", ms = marker_size, capsize=cap_size, ls=ls_err, lw=lw_err, capthick=lw_err, zorder=3)
	axins2.scatter(xN, fidel_deph_true, marker = "x", lw=lw_err, label = "true fidelity (phase free)", zorder=2, color="darkgreen")
	# zoom in on N=2 data for a dephased initial state

	axs[1].add_artist(zoom_axs1)
	axs[2].add_artist(zoom_axs2)

	axs[1].legend(ncol = 4, columnspacing = 0.3, loc = "center",
	  		  handlelength = 1.25, labelspacing = 0.5, framealpha = 0.8, borderpad = 0.3,
	  		  borderaxespad = 0.2, handletextpad = 0.2, fontsize = 8, bbox_to_anchor = (0.5, -0.16), )
	
	plt.tight_layout()
	plt.subplots_adjust(hspace=0.3)
	name = "N_scaling"
	PyPlot.savefig(string(path, name, ".svg"), transparent=true)
	PyPlot.savefig(string(path, name, ".pdf"))
	PyPlot.close(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═b7bf6320-ae60-11ef-06f3-2d4cb0982e2d
# ╠═38280dd2-5a30-46a8-8f08-a3aef31dc081
# ╠═a9fd6f6a-9461-4b4d-8e11-0c64d7e1918d
# ╠═a9bdea14-43a0-4550-bee3-4691141e19e3
# ╠═46ad260d-44ae-4eb3-afb9-2a16e72e2e2f
# ╠═1f8b29b6-03f0-46fa-a455-ca789a4725d0
# ╠═537b1e3e-d4e7-4bd3-872d-83b9fbd3c15f
# ╠═3d043162-1554-4ef1-9709-73bfacc31346
# ╠═e1027453-7355-4ceb-8cbf-6fbcb89e629b
