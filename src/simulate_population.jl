#!/usr/bin/env julia

## simulate_population.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Uses WaveFit to perform simple valley crossing simulations.

#using Revise
include("WaveFit.jl")
using .WaveFit
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
            default=10.0
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
            arg_type =Int64
            default=100
        "--file"
            arg_type=AbstractString
            default = "test"
    end

    parsed_args = parse_args(args, s)

    K = round.(Int64,parsed_args["K"])

    job_name, sigma, delta, s, UL, beta, mu1, mu2, file, num_crossings =
        parsed_args["job_name"], parsed_args["sigma"], parsed_args["delta"],
        parsed_args["s"], parsed_args["UL"], parsed_args["beta"],
        parsed_args["mu1"], parsed_args["mu2"], parsed_args["file"],
        parsed_args["num_crossings"]

    params = parsed_args

    burn_time = K/10

    landscape = Landscape(sigma, delta, s, UL, beta, [0.0, 0.0])
    pop = Population(K,landscape)

    outfile = file

    crossing_times = []

    df = DateFormat("yyyy.mm.dd HH:MM:SS")
    println(Dates.format(now(), df), " | start")
    for i in 1:num_crossings
        pop = Population(K,landscape)
        while get_frequencies(pop)[2] < 0.5
            evolve_multi!(pop)
            if pop.generation == burn_time
                pop.landscape = Landscape(sigma, delta, s, UL, beta, [mu1, mu2])
            end
        end
        push!(crossing_times, pop.generation-burn_time)
        println(Dates.format(now(), df), " | crossing time: ", pop.generation-burn_time)
        flush(stdout)
        file = occursin(r"\.jld2$", outfile) ? outfile : outfile*".jld2"
        @save "output/vc_sims/$file" crossing_times params
    end
end

#main(ARGS)
