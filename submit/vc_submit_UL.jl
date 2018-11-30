#!/usr/bin/env julia

## hoc_submit.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Submit script for simulate_population.jl;
## runs valley crossing sims in batch using SLURM.

using JLD

pars = ["sigma", "K", "mu", "UL"]
sigmalist = collect(logspace(-5,-1,9))
#NGammalist = collect(linspace(.01,.15,8))

#sulist = 5.0*collect(logspace(-4,-2,10))
Klist = [1e5]
mulist = [3e-6]
ULlist = [0.1,1,10,100]
parvals = [sigmalist, Klist, mulist, ULlist]
# take the Cartesian product of all parameter combinations
parsets = collect(Base.product(parvals...))

nsets = length(parsets)

basename = "neut_tunnel_UL"

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
