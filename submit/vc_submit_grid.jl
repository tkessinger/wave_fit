#!/usr/bin/env julia

## vc_submit_grid.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Submit script for simulate_population.jl;
## runs valley crossing sims in batch using SLURM.
## Produces a grid of valley crossing times versus N and δ.

using JLD

pars = ["sigma", "K", "mu"]
#sigmalist = collect(logspace(-5,-1,9))
sigmalist = [1e-2,1e-6]
#NGammalist = collect(linspace(.01,.15,8))

#sulist = 5.0*collect(logspace(-4,-2,10))
Klist = collect(logspace(3,6,7))
mulist = [1e-5,5e-5,1e-4]
slist = [1e-1]
deltalist = collect(logspace(-5,-2,7))
push!(deltalist, 0)
parvals = [sigmalist, Klist, mulist, deltalist]
# take the Cartesian product of all parameter combinations
parsets = collect(Base.product(parvals...))

nsets = length(parsets)

basename = "grid_tunnel"

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
