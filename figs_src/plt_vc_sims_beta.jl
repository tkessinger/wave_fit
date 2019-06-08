#!/usr/bin/env julia

## plt_vc_sims.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Open outputs of simulate_population.jl and plot results.
## This yields plots of crossing time versus Nσ for variable UL.

#ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"

using PyPlot, JLD, Glob, WaveFit, LaTeXStrings, StatsBase

date = "20181219"

    files = glob("output/vc_sims/" * date * "_var_beta_*_*.jld2")
    # glob all the output files

    crossing_times = Dict()

    for file in files
        f = load(file)
        #println(f["params"])
        params = f["parsed_args"]
        sigma, beta, N = params["sigma"], params["beta"], params["K"]
        if (sigma, beta, N) in keys(crossing_times)
            [push!(crossing_times[(sigma, beta, N)], x) for x in f["crossing_times"]]
        else
            crossing_times[(sigma, beta, N)] = f["crossing_times"]
        end
        println("$file, $sigma, $beta, $N, $(length(f["crossing_times"])), $(mean(f["crossing_times"]))")
    end
    sigmalist = sort(unique([x[1] for x in keys(crossing_times)]))
    betalist = sort(unique([x[2] for x in keys(crossing_times)]))
    Nlist = sort(unique([x[3] for x in keys(crossing_times)]))

    tau_grid = zeros(length(sigmalist), length(betalist), length(Nlist))

    for (si, sigma) in enumerate(sigmalist)
        for (bi, beta) in enumerate(betalist)
            for (Ni, N) in enumerate(Nlist)
                if (sigma, beta, N) in keys(crossing_times)
                    tau_grid[si,bi,Ni] = mean(crossing_times[(sigma,beta,N)])
                end
            end
        end
    end

figure()
for (si, sigma) in enumerate(sigmalist)
    for (bi, beta) in enumerate(betalist)
        if abs(beta) == 0.1
        ax = errorbar(Nlist, [mean(x) for x in tau_grid[si,bi,:]],
            yerr = [std(x) for x in tau_grid[si,bi,:]],
            label = latexstring("\$ \\beta = $beta, \\sigma = 10^{$(round(Int64,log10(sigma)))} \$"))
        end
    end
end

ax = axes()
legend(loc=2)
ylabel(L"$\tau$")
xlabel(L"$N$")
xscale("log")
yscale("log")
ylim([1e2,1e7])
tight_layout()
show()
#savefig("neut_tunnel.pdf")
