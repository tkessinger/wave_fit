#!/usr/bin/env julia

## plt_grid_sims.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Make "grid" plots of crossing time versus δ and N.

using PyPlot, JLD, Glob, LaTeXStrings #WaveFit, LaTeXStrings
#gr()

#function load_and_plot(namestring::AbstractString)

date = "20181210"

    files = glob("output/vc_sims/" * date * "_grid_tunnel_*_*.jld")
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
        println("$file, $sigma, $delta, $N, $(length(f["crossing_times"])), $(mean(f["crossing_times"]))")
    end
    sigmalist = sort(unique([x[1] for x in keys(crossing_times)]))
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

    logNlist = log10.(Nlist)
    xticklabels = append!([L"0"],latexstring.(["10^{"*string(round(x*4)/4)*"}" for x in logNlist]))
    logdeltalist = log10.(deltalist[2:end])
    yticklabels = append!([L"0",L"0"],latexstring.(["10^{"*string(round(x*4)/4)*"}" for x in logdeltalist]))


    fig, ax1 = subplots(1,1)
    imshow(max.(0,log10.(tau_grid[1,:,:])),vmin=0,vmax=6)
    cbar = colorbar()
    ax1[:set_xticklabels](xticklabels)
    ax1[:set_yticklabels](yticklabels)
    cbar[:set_label](L"$\log_10 \tau$")
    title(L"$\sigma = 10^{-6}$")

    fig, ax2 = subplots(1,1)
    imshow(max.(0,log10.(tau_grid[2,:,:])),vmin=0,vmax=6)
    cbar = colorbar()
    ax2[:set_xticklabels](xticklabels)
    ax2[:set_yticklabels](yticklabels)
    cbar[:set_label](L"$\log_10 \tau$")
    title(L"$\sigma = 10^{-2}$")


    #println(param_combs)
#end
#sigmalist, UL_list, nt_results = load_and_plot("neut_tunnel")
#UL_vals = unique(UL_list)
#display(
#load_and_plot("grid_tunnel")
