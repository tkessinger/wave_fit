#!/usr/bin/env julia

## simulate_population.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Uses WaveFit to perform simple valley crossing simulations.

#using Revise
using WaveFit
using Distributions
using Stats
using ArgParse, JLD


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

    parsed_args = parse_args(s)

    K = floor.(Int64,parsed_args["K"])

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

    for i in 1:num_crossings
        pop = Population(K,landscape)
        while get_frequencies(pop)[2] < 0.5
            evolve!(pop)
            if pop.generation == burn_time
                pop.landscape = Landscape(sigma, delta, s, UL, beta, [mu1, mu2])
            end
        end
        push!(crossing_times, pop.generation-burn_time)
        file = ismatch(r"\.jld", outfile) ? outfile : outfile*".jld"
        save("output/vc_sims/$file", "crossing_times", crossing_times, "params", params)
    end
end

main(ARGS)
