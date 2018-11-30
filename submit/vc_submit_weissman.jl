#!/usr/bin/env julia

## vc_submit_weissman.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Submit script for simulate_population.jl;
## runs valley crossing sims in batch using SLURM.
## Duplicates main figures from Weissman et al. (2009)

using JLD

fig_to_make = 4

function weissman_fig(fig_to_make)
        slist = [0.1]
        sigmalist = [1e-2,1e-8]
        if fig_to_make == 1
                Klist = collect(logspace(2,6,9))
                mu1list = [1e-5]
                mu2list = [1e-4]
                deltalist = [2e-4]
        elseif fig_to_make == 2
                Klist = collect(logspace(2,6,9))
                mu1list = [1e-5]
                mu2list = [1e-6]
                deltalist = [7e-3]
        elseif fig_to_make == 3
                Klist = [1e5]
                mu1list = [1e-6]
                mu2list = [4e-5]
                deltalist = [-1e-3,-5e-3,-1e-3,0,1e-3,5e-3,1e-2]
        elseif fig_to_make == 4
                Klist = [1e5]
                mu1list = [1e-6]
                mu2list = collect(logspace(-6,-2,9))
                deltalist = [3e-3]
        end

        # take the Cartesian product of all parameter combinations
        pars = ["sigma", "K", "mu1", "mu2", "s", "delta"]
        parvals = [sigmalist, Klist, mu1list, mu2list, slist, deltalist]
        parsets = collect(Base.product(parvals...))

        nsets = length(parsets)

        basename = "weissman_fig_"*string(fig_to_make)

        for p in 1:nsets

                simstr = "src/simulate_population.jl"

                # create filname base from basename and run number
                numstr = lpad(string(p), length(string(nsets)), "0")
                filebase = basename * "_" * numstr

                # add parameter values as command line options
                simstr *= " "
                simstr *= join(map( (x) -> "--" * x[1] * " " * string(x[2]), zip(pars, parsets[p])), " ")
                simstr *= " --file " * filebase * ".jld"

                # run the scripts via sbatch
                # the "split" wizardry is needed due to the fine distinction between Cmd and String types
                print("sbatch $(split(simstr))\n")
                run(`sbatch $(split(simstr))`)
        end
end

weissman_fig(fig_to_make)
