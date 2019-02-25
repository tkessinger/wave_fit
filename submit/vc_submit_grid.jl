#!/usr/bin/env julia

## vc_submit_grid.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Submit script for simulate_population.jl;
## runs valley crossing sims in batch using SLURM.
## Produces a grid of valley crossing times versus N and δ.

using JLD, Dates

date = Dates.format(Dates.Date(Dates.now()), "yyyymmdd")

pars = ["sigma", "K", "mu1", "mu2", "delta", "num_crossings"]
#sigmalist = collect(logspace(-5,-1,9))
sigmalist = [1e-2,1e-6]
#NGammalist = collect(linspace(.01,.15,8))

jobids = []
cmd_strings = []

#sulist = 5.0*collect(logspace(-4,-2,10))
Klist = round.(collect(logspace(3,6,7)),-1)
#don't ask me why the -1 is necessary
mu1list = [1e-5]
mu2list = [1e-5]
slist = [1e-1]
deltalist = collect(logspace(-5,-2,7))
push!(deltalist, 0)
num_crossings_list = [10]
parvals = [sigmalist, Klist, mu1list, mu2list, deltalist, num_crossings_list]
# take the Cartesian product of all parameter combinations
parsets = collect(Base.product(parvals...))

nsets = length(parsets)

runs_per_param_comb = 10

basename = date*"_grid_tunnel"

for p in 1:nsets
    for runno in 1:runs_per_param_comb

        simstr = "src/simulate_population.jl"

        # create filname base from basename and run number
        numstr = lpad(string(p), length(string(nsets)), "0")
        runstr = lpad(string(runno), 2, "0")
        filebase = basename * "_" * numstr * "_" * runstr

        # add parameter values as command line options
        simstr *= " "
        simstr *= join(map( (x) -> "--" * x[1] * " " * string(x[2]), zip(pars, parsets[p])), " ")
        simstr *= " --file " * filebase * ".jld"

        # run the scripts via sbatch
        # the "split" wizardry is needed due to the fine distinction between Cmd and String types
        println(`sbatch $(split(simstr))`)
        jobout = String(read(`sbatch $(split(simstr))`))

        #jobout = String(read(`sbatch $simstr`))

        push!(cmd_strings, `sbatch $(split(simstr))`)
        push!(jobids, match(r"[0-9]+", jobout).match)
    end
end

jid_file = "$(date)_grid_sims_jobids.dat"
cmd_file = "$(date)_grid_sims_commands.dat"

writedlm(jid_file, [[(jobids[p], p, runno) for runno in runs_per_param_comb] for p in 1:nsets])
println("\njob IDs written to file $jid_file\n")
writedlm(cmd_file, [(jobids[x], cmd_strings[x]) for x in length(jobids)])
println("submitted commands written to file $cmd_file\n")
