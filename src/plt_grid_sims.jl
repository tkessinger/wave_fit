#!/usr/bin/env julia

## plt_grid_sims.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Make "grid" plots of crossing time versus δ and N.

using PyPlot, JLD, Glob, WaveFit, LaTeXStrings
#gr()

#function load_and_plot(namestring::AbstractString)

    files = glob("output/vc_sims/grid_tunnel_*.jld")
    # glob all the output files

    crossing_times = Dict()

    for file in files
        f = load(file)
        #println(f["params"])
        params = f["params"]
        sigma, delta, N = params["sigma"], params["delta"], params["K"]
        if (sigma, delta, N) in keys(crossing_times)
            [push!(crossing_times[(sigma, delta, N)], x) for x in f["crossing_times"]]
        else
            crossing_times[(sigma, delta, N)] = f["crossing_times"]
        end
        println("$file, $sigma, $delta, $N, $(mean(f["crossing_times"]))")
    end
    sigmalist = unique([x[1] for x in keys(crossing_times)])
    deltalist = sort(unique([x[2] for x in keys(crossing_times)]))
    Nlist = sort(unique([x[3] for x in keys(crossing_times)]))

    tau_grid = zeros(length(sigmalist), length(deltalist), length(Nlist))

    for (si, sigma) in enumerate(sigmalist)
        for (di, delta) in enumerate(deltalist)
            for (Ni, N) in enumerate(Nlist)
                if (sigma, delta, N) in keys(crossing_times)
                    tau_grid[si,di,Ni] = mean(crossing_times[(sigma,delta,N)])
                end
            end
        end
    end

    imshow(tau_grid)

    #println(param_combs)
#end
#sigmalist, UL_list, nt_results = load_and_plot("neut_tunnel")
#UL_vals = unique(UL_list)
#display(
#load_and_plot("grid_tunnel")
