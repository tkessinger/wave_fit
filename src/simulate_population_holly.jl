#!/usr/bin/julia

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=1-00:00:00

# Set name of job shown in squeue
#SBATCH --job-name sample_gaussian.jl

# Request CPU resources
#SBATCH --qos=max40
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=2
#SBATCH --account=8e300-jva238-lowprio
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Uses WaveFit to perform simple valley crossing simulations.

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
            default=1e5
        "--sigma"
            arg_type=Float64
            default=1e-5
        "--delta"
            arg_type=Float64
            default=1e-3
        "--s"
            arg_type=Float64
            default=1e-2
        "--UL"
            arg_type=Float64
            default=1.0
        "--beta"
            arg_type=Float64
            default=0.01
        "--mu1"
            arg_type=Float64
            default=1e-5
        "--mu2"
            arg_type=Float64
            default=1e-5
        "--num_crossings"
            arg_type=Int64
            default=100
        "--burn_factor"
            arg_type=Float64
            default=0.1
        "--file"
            arg_type=AbstractString
            default = "test"
    end

    seed = Base.Random.GLOBAL_RNG.seed

    parsed_args = parse_args(args, s)

    job_name, K, sigma, delta, s, UL, beta, mu1, mu2, outfile, num_crossings, burn_factor =
        parsed_args["job_name"], parsed_args["K"], parsed_args["sigma"],
        parsed_args["delta"], parsed_args["s"], parsed_args["UL"], parsed_args["beta"],
        parsed_args["mu1"], parsed_args["mu2"], parsed_args["file"],
        parsed_args["num_crossings"], parsed_args["burn_factor"]

    println("-----------")
    foreach(p -> println(p[1], ": ", p[2]), parsed_args)
    println("random seed: $seed")
    println("-----------")

    flush(STDOUT)

    K = round(Int64, K)
    burn_time = round(Int64, K * burn_factor)
    landscape = Landscape(sigma, delta, s, beta, UL, [0.0, 0.0])
    pop = Population(K,landscape)

    crossing_times = []

    df = DateFormat("yyyy.mm.dd HH:MM:SS")
    println(Dates.format(now(), df), " | start")
    for i in 1:num_crossings
        pop = Population(K,landscape)
        while get_frequencies(pop)[2] < 0.5
            evolve_multi!(pop)
            if pop.generation == burn_time
                pop.landscape = Landscape(sigma, delta, s, beta, UL, [mu1, mu2])
            end
        end
        push!(crossing_times, pop.generation-burn_time)
        println(Dates.format(now(), df), " | crossing time: ", pop.generation-burn_time)
        flush(STDOUT)
        file = ismatch(r"\.jld2", outfile) ? outfile : outfile*".jld2"
        @save "/scratch/tke235/output/vc_sims/$file" crossing_times parsed_args seed
    end
end

main(ARGS)
