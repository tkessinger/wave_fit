#!/usr/bin/env julia

## simulate_population.jl
##
## Use WaveFit to perform simple valley crossing simulations.
## Run simulations on via multicore processesing.

using WaveFit
using Distributions
using ArgParse
using Dates
using Distributed
using DataFrames
using CSV
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

function main(args)

    s = ArgParseSettings(description = "Run WaveFit simulation across multiple cores")
    @add_arg_table s begin
        "--ncpus"
            arg_type=Int64
            default=max(round(Int,Sys.CPU_THREADS/2), 1)
        "--input"
            default=nothing
        "--output"
            arg_type=String
            default="test"
    end
    parsed_args = parse_args(args, s)

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
    pars = read_parameters(defpars, parsed_args["input"])

    # take the Cartesian product of all parameter combinations
    parsets = collect(Base.product(values(pars)...))
    nsets = length(parsets)

    # setup workers assuming directory is manually added to LOAD_PATH
    addprocs(min(parsed_args["ncpus"], round(Int64, Sys.CPU_THREADS/2)))
    extradir = filter((p)->match(r"/", p) !== nothing, LOAD_PATH)[1]
    @everywhere workers() push!(LOAD_PATH, $extradir)
    @everywhere workers() eval(:(using Random))
    @everywhere workers() eval(:(using WaveFit))
    @everywhere workers() eval(:(using Dates))

    inputs  = RemoteChannel(()->Channel{Dict}(2*nsets*maximum(pars["num_crossings"])))
    results = RemoteChannel(()->Channel{Dict}(2*nsets*maximum(pars["num_crossings"])))

    @everywhere function run_worker(inputs, results)
        # save crossing number and random seed
        seed = Dict(zip(["seed1", "seed2", "seed3", "seed4"], Random.GLOBAL_RNG.seed))

        while true
            pard = take!(inputs)
            pard = merge(pard, seed)
            burn_time = round(Int64, pard["K"] * pard["burn_factor"])

            println("--- running ", pard["nrun"], " --- ")
            flush(stdout)

            # burn in
            start = now()
            landscape = Landscape(pard["sigma"],
                                  pard["delta"],
                                  pard["s"],
                                  pard["beta"],
                                  pard["UL"],
                                  [0.0, 0.0])

            pop = Population(pard["K"], landscape)
            for i in 1:burn_time
                evolve_multi!(pop)
            end

            # initialize mutations and run until valley crossing
            pop.landscape = Landscape(pard["sigma"],
                                      pard["delta"],
                                      pard["s"],
                                      pard["beta"],
                                      pard["UL"],
                                      [pard["mu1"], pard["mu2"]])
            while get_frequencies(pop)[2] < 0.5
                evolve_multi!(pop)
            end

            # save crossing time
            pard["crossing_time"] = pop.generation-burn_time

            # output elapsed time
            stop = now()
            print("--- ran ", pard["nrun"], " --- elapsed time: ",
                Dates.canonicalize(Dates.CompoundPeriod(round(stop-start, Dates.Second(1)))), " ")
            foreach(k -> print(k, ": ", pard[k], ", "),
                sort(collect(keys(filter(p->p.first ∉ ["nrun"], pard)))))
            println()
            flush(stdout)

            # return data to master process
            put!(results, pard)
        end
    end

    # load parameter sets into inputs channel
    nruns = 0
    for parset in parsets
        pard = Dict(zip(keys(pars), parset))
        print("--- queueing --- ")
        foreach(k -> print(k, ": ", pard[k], ", "), sort(collect(keys(pard))))
        println()
        flush(stdout)
        for rep in 1:pard["num_crossings"]
            nruns += 1
            rpard = copy(pard)
            rpard["rep"] = rep
            rpard["nrun"] = nruns
            put!(inputs, rpard)
        end
    end

    # start workers running on parameter sets in inputs
    for w in workers() # start tasks on the workers to process requests in parallel
        remote_do(run_worker, w, inputs, results)
    end

    # create output file name and data table
    output = parsed_args["output"]
    file = occursin(r"\.csv$", output) ? output : output*".csv"
    cols = push!(sort(collect(keys(pars))),
                 ["rep", "crossing_time", "seed1", "seed2", "seed3", "seed4"]...)
    dat = DataFrame(Dict([(c, Any[]) for c in cols]))

    # grab results and output to CSV
    for sim in 1:nruns
        # get results from parallel jobs
        flush(stdout)
        resd = take!(results)
        nrun = pop!(resd, "nrun")

        # add to table (must convert dict keys to symbols) and save
        push!(dat, Dict([(Symbol(k), resd[k]) for k in keys(resd)]))
        CSV.write(file, dat)
    end
end

main(ARGS)
