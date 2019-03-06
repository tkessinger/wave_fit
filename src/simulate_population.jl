#!/usr/bin/env julia

## simulate_population.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Uses WaveFit to perform simple valley crossing simulations.

#using Revise
#include("WaveFit.jl")
using WaveFit
using Distributions
using ArgParse, JLD2
using Dates

function main(args)

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
            arg_type=Float64
            default=100
        "--burn_factor"
            arg_type=Float64
            default=0.25
        "--file"
            arg_type=AbstractString
            default = "test"
    end

    parsed_args = parse_args(args, s)

    job_name, K, sigma, delta, s, UL, beta, mu1, mu2, outfile, num_crossings, burn_factor =
        parsed_args["job_name"], parsed_args["K"], parsed_args["sigma"],
        parsed_args["delta"], parsed_args["s"], parsed_args["UL"], parsed_args["beta"],
        parsed_args["mu1"], parsed_args["mu2"], parsed_args["file"],
        parsed_args["num_crossings"], parsed_args["burn_factor"]

    println("-----------")
    foreach(p -> println(p[1], ": ", p[2]), parsed_args)
    println("-----------")

    K = round(Int64, K)
    burn_time = round(Int64, K * burn_factor)
    num_crossings = round(Int64, num_crossings)
    landscape = Landscape(sigma, delta, s, beta, UL, [0.0, 0.0])

    crossing_times = []

    df = DateFormat("yyyy.mm.dd HH:MM:SS")
    println(Dates.format(now(), df), " | start")
    for i in 1:num_crossings
        pop = Population(K, landscape)

        # burn in
        for i in 1:burn_time
            evolve_multi!(pop)
        end

        # initialize mutations and run until valley crossing
        pop.landscape = Landscape(sigma, delta, s, beta, UL, [mu1, mu2])
        while get_frequencies(pop)[2] < 0.5
            evolve_multi!(pop)
        end

        push!(crossing_times, pop.generation-burn_time)
        println(Dates.format(now(), df), " | crossing time: ", pop.generation-burn_time)
        flush(stdout)
        file = occursin(r"\.jld2$", outfile) ? outfile : outfile*".jld2"
        @save "$file" crossing_times parsed_args
    end
end

main(ARGS)
