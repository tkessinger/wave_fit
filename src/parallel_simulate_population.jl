#!/usr/bin/env julia

## simulate_population.jl
##
## Use WaveFit to perform simple valley crossing simulations.
## Run simulations on via multicore processesing.

using WaveFit
using Distributions
using ArgParse, JLD2
using Dates
import JSON

# read parameters from JSON file
function read_parameters(defpars, inputfile=nothing)
    pars = copy(defpars)

    # read JSON file
    if inputfile != nothing
        inpars = JSON.parsefile(inputfile)
    else
        inpars = Dict()
    end

    for parkey in keys(defpars)
        if "type" in keys(pars[parkey]) &&
            isprimitivetype(pars[parkey]["type"])
            T = pars[parkey]["type"]
        else
            # Default type is Float
            T = Float64
        end
        if T <: Int
            convertf = (val) -> round(T, val)
        else
            convertf = (val) -> convert(T, val)
        end

        # use defpars for list of usabe parameters in JSON
        if parkey in keys(inpars)
            if "value" in keys(inpars[parkey])
                val = inpars[parkey]["value"]
            elseif "range" in keys(inpars[parkey])
                valr = inpars[parkey]["range"]
                if "log" in keys(valr)
                    b = valr["log"]
                    rf = (r) -> b .^ r
                    pop!(valr, "log")
                else
                    rf = (r) -> r
                end
                start = pop!(valr, "start")
                rkws = Dict(zip(Symbol.(keys(valr)), values(valr)))
                val = rf(range(start; rkws...))
            end
        else
            val = pars[parkey]["value"]
        end

        if !isstructtype(typeof(val))
            pars[parkey] = [convertf(val)]
        else
            pars[parkey] = convertf.(val)
        end
    end

    return pars
end

pars = read_parameters(defpars, "test_parameter_values.json")

function main(args)

    defpars = Dict{String,Any}([
        "K"     => Dict("value" => 10^5, "type" => Int64),
        "sigma" => Dict("value" => 1e-5, "type" => Float64),
        "delta" => Dict("value" => 1e-3, "type" => Float64),
        "s"     => Dict("value" => 1e-2, "type" => Float64),
        "UL"    => Dict("value" => 1.0, "type" => Float64),
        "beta"  => Dict("value" => 0.01, "type" => Float64),
        "mu1"   => Dict("value" => 1e-5, "type" => Float64),
        "mu2"   => Dict("value" => 1e-5, "type" => Float64),
        "num_crossings" => Dict("value" => 100, "type" => Int64),
        "burn_factor"   => Dict("value" => 0.25, "type" => Float64)
    ])

    s = ArgParseSettings(description = "Run WaveFit simulation across multiple cores")

    @add_arg_table s begin
        "--job_name"
            arg_type=String
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
            default=100.0
        "--burn_factor"
            arg_type=Float64
            default=0.25
        "--input"
            default=nothing
        "--output"
            arg_type=String
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
