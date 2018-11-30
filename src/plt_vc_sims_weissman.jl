#!/usr/bin/env julia

## plt_vc_sims_weissman.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Reproduces main crossing time plots from Weissman et al. (2009).


using PyPlot, JLD, Glob, WaveFit, LaTeXStrings

function load_and_plot(fig_to_make::Int64)
    if fig_to_make==1 || fig_to_make==2
        indep_var = "K"
        indep_var_string = L"N"
    elseif fig_to_make==3
        indep_var = "delta"
        indep_var_string = L"\delta"
    elseif fig_to_make==4
        indep_var = "mu2"
        indep_var_string = L"\mu_1"
    end
    sigmalist = [1e-8,1e-2]
    low_sigma_results = []
    high_sigma_results = []

    namestring = "weissman_fig_"*string(fig_to_make)
    files = glob("output/vc_sims/" * namestring * "*.jld")
    # glob all the output files

    for (fi, file) in enumerate(files)
        f = load(file)
        file_results = convert(Array{Int64,1},f["crossing_times"])
        indep_var_vals = f["params"][indep_var]
        println("$fig_to_make, $indep_var_vals, $(length(file_results))")
        if f["params"]["sigma"] == sigmalist[1]
            push!(low_sigma_results, [indep_var_vals, file_results])
        elseif f["params"]["sigma"] == sigmalist[2]
            push!(high_sigma_results, [indep_var_vals, file_results])
        end
    end
    #println(low_sigma_results)
    #println(high_sigma_results)

    fig = figure()
    ax1 = errorbar(
        [x[1] for x in low_sigma_results],
        [mean(x[2]) for x in low_sigma_results],
        yerr = [std(x[2]) for x in low_sigma_results],
        label = L"\sigma = 10^{-8}"
    )
    ax2 = errorbar(
        [x[1] for x in high_sigma_results],
        [mean(x[2]) for x in high_sigma_results],
        yerr = [std(x[2]) for x in high_sigma_results],
        label = L"\sigma = 10^{-2}"
    )
    ylabel(L"\tau")
    xlabel(indep_var_string)
    title("figure \($(collect('a':'z')[fig_to_make])\)")
    if (fig_to_make != 3)
        xscale("log")
        yscale("log")
    end
    legend(loc=2)
    savefig("julia_weissman_$fig_to_make.pdf")
#    return low_sigma_results, high_sigma_results

#    return (sigmalist, UL_list, all_results)
end

#(low_sigma_results, high_sigma_results) = load_and_plot(1)
load_and_plot(1)
load_and_plot(2)
load_and_plot(3)
load_and_plot(4)

# sigmalist, UL_list, nt_results = load_and_plot("neut_tunnel")
# UL_vals = unique(UL_list)
# #display(
# fig = figure()
#
# for (ii, i) in enumerate(UL_vals)
#     ax = errorbar(
#         sigmalist,
#         [mean(x) for x in nt_results[((ii-1)*length(sigmalist)+1):ii*length(sigmalist)]],
#         yerr = [1.0*std(x) for x in nt_results[((ii-1)*length(sigmalist)+1):ii*length(sigmalist)]],
#         label = "UL = $i"
#         )
# end
# ax = axes()
# legend(loc=2)
# ylabel(L"$\tau$")
# xlabel(L"$\sigma$")
# xscale("log")
# yscale("log")
# title("neutral tunneling")
# ylim([1e4,1e7])
# #savefig("neut_tunnel.pdf")
