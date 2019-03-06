#!/usr/bin/env julia

## plt_vc_sims_weissman.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Reproduces main crossing time plots from Weissman et al. (2009).


using Statistics, StatsBase, PyPlot, FileIO, JLD2, Glob, WaveFit, LaTeXStrings

function load_and_plot(fig_to_make::Int64, equalsize=false)
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
    sigmalist = [0.5,1e-1,0.05,1e-2,1e-8]
    results = Dict()

    namestring = "weissman_fig_"*string(fig_to_make)
    files = glob(namestring * "*.jld2")
    # glob all the output files

    for (fi, file) in enumerate(files)
        ctimes, pargs = load(file, "crossing_times", "parsed_args")
        file_results = convert(Array{Int64,1}, ctimes)
        indep_var_vals = pargs[indep_var]
        # println("$fig_to_make, $indep_var_vals, $(length(file_results))")
        if !(pargs["sigma"] in sigmalist)
            println("$(pargs[sigma]) not in sigmalist")
            continue
        end
        if pargs["sigma"] in keys(results)
            push!(results[pargs["sigma"]], [indep_var_vals, file_results])
        else
            results[pargs["sigma"]] = []
        end
    end
    minsize = minimum([minimum([length(x[2]) for x in results[sigma]]) for sigma in keys(results)])
    if equalsize == 1
        println("Number of samples: $minsize")
    end

    fig = figure()
    ax = []
    for sigma in sort(collect(keys(results)), rev=true)
        len = length(results[sigma])
        if equalsize <= 0
            samps = (:)
        elseif equalsize == 1
            samps = 1:minsize
        elseif equalsize < minsize
            samps = 1:equalsize
        else
            error("some sample sizes too small for " * string(equalsize) * " samples")
        end

        push!(ax,
            errorbar(
                     [x[1] for x in results[sigma]],
                     [median(x[2][samps]) for x in results[sigma]],
                     yerr = vcat([quantile(x[2][samps], 0.25) for x in results[sigma]]',
                                 [quantile(x[2][samps], 0.75) for x in results[sigma]]'),
                     fmt = "o",
                     alpha = 0.75,
                     label = L"$\sigma$ = " * string(sigma),
                     )
            )
    end

    ylabel(L"\tau")
    xlabel(indep_var_string)
    title("figure ($(collect('a':'z')[fig_to_make]))")
    if (fig_to_make != 3)
        xscale("log")
        yscale("log")
    end
    legend(loc=3)
    savefig("julia_weissman_$fig_to_make.pdf")
#    return low_sigma_results, high_sigma_results

#    return (sigmalist, UL_list, all_results)
end

#(low_sigma_results, high_sigma_results) = load_and_plot(1)
load_and_plot(1, false)
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
