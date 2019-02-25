#!/usr/bin/env julia

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=1-00:00:00

# Set name of job shown in squeue
#SBATCH --job-name sample_gaussian.jl

# Request CPU resources
#SBATCH --qos=max40
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

## sample_gaussian.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Uses WaveFit and records the number and size of fitness classes.

#using Revise
using WaveFit
using Distributions
using Stats, Dates
using ArgParse, JLD2


function main(args)

    @show args

    s = ArgParseSettings(description = "Example 1 for argparse.jl: minimal usage.")

    @add_arg_table s begin
        "--job_name"
            arg_type=AbstractString
            default = "test"
        "--K"
            arg_type=Float64
            default=1e4
        "--sigma"
            arg_type=Float64
            default=1e-5
        "--UL"
            arg_type=Float64
            default=10.0
        "--beta"
            arg_type=Float64
            default=0.01
        "--num_samples"
            arg_type=Int64
            default=10
        "--sampling_interval"
            arg_type=Float64
            default=0.01
        "--file"
            arg_type=AbstractString
            default = "test"
    end

    seed = Base.Random.GLOBAL_RNG.seed

    parsed_args = parse_args(args, s)

    job_name, K, sigma, UL, beta, outfile, num_samples, sampling_interval =
        parsed_args["job_name"], parsed_args["K"], parsed_args["sigma"],
        parsed_args["UL"], parsed_args["beta"],
        parsed_args["file"],
        parsed_args["num_samples"],
        parsed_args["sampling_interval"]

    println("-----------")
    foreach(p -> println(p[1], ": ", p[2]), parsed_args)
    println("-----------")

    flush(STDOUT)

    K = round(Int64, K)
    sampling_time = round(Int64, K * sampling_interval)
    landscape = Landscape(sigma, 0.0, 0.0, beta, UL, [0.0, 0.0])
    pop = Population(K,landscape)

    pop_samples = []
    samples = 0

    clones = []

    df = DateFormat("yyyy.mm.dd HH:MM:SS")
    println(Dates.format(now(), df), " | start")
    pop = Population(K,landscape)
    while samples < num_samples
        evolve_multi!(pop)
        if (pop.generation%sampling_time == 0)
            clones = [(x.bg_mutations, x.n) for x in collect(values(pop.classes))]
            samples += 1
            push!(pop_samples, clones)
            println(Dates.format(now(), df), " | sampled at: ", pop.generation, ", $(length(clones)) clones")
            flush(STDOUT)
            file = ismatch(r"\.jld2", outfile) ? outfile : outfile*".jld2"
            println("$(length(pop_samples))")
            @save "output/vc_sims/$file" pop_samples parsed_args seed
        end
    end
    return pop_samples
end

ps = main(ARGS)
