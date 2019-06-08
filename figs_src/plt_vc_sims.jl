#!/usr/bin/env julia

## plt_vc_sims.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Open outputs of simulate_population.jl and plot results.

ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"

using Plots, JLD, Glob, WaveFit, LaTeXStrings
gr()

function load_and_plot(namestring::AbstractString)
    if namestring == "seq_fix"
        K = floor.(Int64, 1e4)
    elseif namestring == "neut_tunnel"
        K = floor.(Int64, 1e5)
    end
    files = glob("output/vc_sims/" * namestring * "_*.jld")
    # glob all the output files

    # keys for this dict will be K value and s_u value
    sigmalist = collect(logspace(-5,-1,9))

    all_results = Array{Int64,1}[]
    # open files and fill dict
    for (fi, file) in enumerate(files)
        #println("$file")
        f = load(file)
        results = convert(Array{Int64,1},f["crossing_times"])
        results -= K/10
        push!(all_results, results)
    end

    #fig = figure()
    #ax = fig[:add_subplot](1,1,1)
    #display(plot(sigmalist, [mean(x) for x in all_results], xaxis=:log, yaxis=:log))

    return all_results
end
sigmalist = collect(logspace(-5,-1,9))
seq_results = load_and_plot("seq_fix")
display(
plot(
    sigmalist,
    [mean(x) for x in seq_results],
    yerror = [1.0*std(x) for x in seq_results],
    xaxis=:log, yaxis=:log,
    ylims=(1e3,1e6),
    ylabel=L"\tau", xlabel=L"\sigma",
    lw = 3,
    legend = false,
    title = "Sequential fixation",
    show = true
    )
)
savefig("seq_fix.pdf")
#display(plot(sigmalist, [log10(mean(x)) for x in seq_results], xaxis=:log, yerror = [log10(std(x)) for x in seq_results]))
nt_results = load_and_plot("neut_tunnel")
display(
plot(
    sigmalist,
    [mean(x) for x in nt_results],
    yerror = [1.0*std(x) for x in nt_results],
    xaxis=:log, yaxis=:log,
    ylims=(1e3,1e7),
    ylabel=L"\tau", xlabel=L"\sigma",
    lw = 3,
    legend = false,
    title = "Stochastic tunneling",
    show = true
    )
)
savefig("neut_tunnel.pdf")
