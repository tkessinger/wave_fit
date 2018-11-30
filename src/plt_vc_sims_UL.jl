#!/usr/bin/env julia

## plt_vc_sims.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Open outputs of simulate_population.jl and plot results.
## This yields plots of crossing time versus Nσ for variable UL.

#ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"

using PyPlot, JLD, Glob, WaveFit, LaTeXStrings
gr()

function load_and_plot(namestring::AbstractString)
    if namestring == "seq_fix"
        K = floor.(Int64, 1e4)
    elseif namestring == "neut_tunnel"
        K = floor.(Int64, 1e5)
    end
    files = glob("output/vc_sims/" * namestring * "_UL_*.jld")
    # glob all the output files

    # keys for this dict will be K value and s_u value
    sigmalist = collect(logspace(-5,-1,9))

    all_results = Array{Int64,1}[]
    UL_list = []
    # open files and fill dict
    for (fi, file) in enumerate(files)
        #println("$file")
        f = load(file)
        results = convert(Array{Int64,1},f["crossing_times"])
        UL = f["params"]["UL"]
        #results -= K/10
        push!(all_results, results)
        push!(UL_list, UL)
    end

    #fig = figure()
    #ax = fig[:add_subplot](1,1,1)
    #display(plot(sigmalist, [mean(x) for x in all_results], xaxis=:log, yaxis=:log))

    return (sigmalist, UL_list, all_results)
end
sigmalist, UL_list, nt_results = load_and_plot("neut_tunnel")
UL_vals = unique(UL_list)
#display(
fig = figure()

for (ii, i) in enumerate(UL_vals)
    ax = errorbar(
        sigmalist,
        [mean(x) for x in nt_results[((ii-1)*length(sigmalist)+1):ii*length(sigmalist)]],
        yerr = [1.0*std(x) for x in nt_results[((ii-1)*length(sigmalist)+1):ii*length(sigmalist)]],
        label = "UL = $i"
        )
end
ax = axes()
legend(loc=2)
ylabel(L"$\tau$")
xlabel(L"$\sigma$")
xscale("log")
yscale("log")
title("neutral tunneling")
ylim([1e4,1e7])
#savefig("neut_tunnel.pdf")
