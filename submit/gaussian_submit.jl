#!/usr/bin/env julia

## gaussian_submit.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Submit script for sample_gaussian.jl;
## generates test population with WaveFit and samples
## their fitness distributions.

using JLD2, Dates

date = Dates.format(Dates.Date(Dates.now()), "yyyymmdd")

pars = ["sigma", "K", "UL", "num_samples"]
#sigmalist = collect(logspace(-5,-1,9))
sigmalist = [5e-2,5e-4,5e-6]
#NGammalist = collect(linspace(.01,.15,8))

jobids = []
cmd_strings = []

#sulist = 5.0*collect(logspace(-4,-2,10))
Klist = round.(collect(logspace(3,6,4)),-1)
#don't ask me why the -1 is necessary
UL_list = [0.1,1.0,10.0]
num_samples_list = [200]
parvals = Any[sigmalist, Klist, UL_list, num_samples_list]
# take the Cartesian product of all parameter combinations
# the Any is needed to keep num_crossings from being forced into a Float64
parsets = collect(Base.product(parvals...))

nsets = length(parsets)

runs_per_param_comb = 3

basename = date*"_sample_gaussian"

for p in 1:nsets
    for runno in 1:runs_per_param_comb

        simstr = "src/sample_gaussian.jl"

        # create filname base from basename and run number
        numstr = lpad(string(p), length(string(nsets)), "0")
        runstr = lpad(string(runno), 2, "0")
        filebase = basename * "_" * numstr * "_" * runstr

        # add parameter values as command line options
        simstr *= " "
        simstr *= join(map( (x) -> "--" * x[1] * " " * string(x[2]), zip(pars, parsets[p])), " ")
        simstr *= " --file " * filebase * ".jld2"

        # run the scripts via sbatch
        # the "split" wizardry is needed due to the fine distinction between Cmd and String types
        println(`sbatch $(split(simstr))`)
        jobout = String(read(`sbatch $(split(simstr))`))

        #jobout = String(read(`sbatch $simstr`))

        push!(cmd_strings, `sbatch $(split(simstr))`)
        push!(jobids, match(r"[0-9]+", jobout).match)
    end
end

jid_file = "$(date)_gaussian_sims_jobids.dat"
cmd_file = "$(date)_gaussian_sims_commands.dat"

writedlm(jid_file, [[(jobids[(p-1)*runs_per_param_comb+runno], p, runno) for runno in 1:runs_per_param_comb] for p in 1:nsets])
println("\njob IDs written to file $jid_file\n")
writedlm(cmd_file, [(jobids[x], cmd_strings[x]) for x in 1:length(jobids)])
println("submitted commands written to file $cmd_file\n")
