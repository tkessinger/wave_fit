#!/usr/bin/env julia

## plt_vc_sims.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Open outputs of simulate_population.jl and plot results.
## This yields plots of crossing time versus Nσ for variable UL.

#ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"

using PyPlot, JLD, Glob, LaTeXStrings, Statistics

date = "20190109"

    files = glob("output/vc_sims/" * date * "_var_UL_*_*.jld2")
    # glob all the output files

    crossing_times = Dict()

    for file in files
        f = load(file)
        #println(f["params"])
        params = f["parsed_args"]
        sigma, UL, N = params["sigma"], params["UL"], params["K"]
        if (sigma, UL, N) in keys(crossing_times)
            [push!(crossing_times[(sigma, UL, N)], x) for x in f["crossing_times"]]
        else
            crossing_times[(sigma, UL, N)] = f["crossing_times"]
        end
        println("$file, $sigma, $UL, $N, $(length(f["crossing_times"])), $(mean(f["crossing_times"]))")
    end
    sigmalist = sort(unique([x[1] for x in keys(crossing_times)]))
    UL_list = sort(unique([x[2] for x in keys(crossing_times)]))
    Nlist = sort(unique([x[3] for x in keys(crossing_times)]))

    tau_grid = zeros(length(sigmalist), length(UL_list), length(Nlist))

    for (si, sigma) in enumerate(sigmalist)
        for (ui, UL) in enumerate(UL_list)
            for (Ni, N) in enumerate(Nlist)
                if (sigma, UL, N) in keys(crossing_times)
                    tau_grid[si,ui,Ni] = mean(crossing_times[(sigma,UL,N)])
                end
            end
        end
    end

for (si, sigma) in enumerate(sigmalist)
    figure()
    for (ui, UL) in enumerate(UL_list)
    ax = errorbar(Nlist, [mean(x) for x in tau_grid[si,ui,:]],
        yerr = [std(x) for x in tau_grid[si,ui,:]],
        label = L"$UL = $" * string(UL))
    end
    legend(loc=2)
    title(L"$\sigma = 10^{" * string(round(Int64,log10(sigma))) * L"}\$")
    ylabel(L"$\tau$")
    xlabel(L"$N$")
    xscale("log")
    yscale("log")
    ylim([1e2,1e7])
    show()
end
