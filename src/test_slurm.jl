#!/usr/bin/env julia

using Distributed
using ClusterManagers

function main(args)

    addprocs(SlurmManager(args[1]))

    hosts = []
	pids = []
	for i in workers()
		host, pid = fetch(@spawnat i (gethostname(), getpid()))
		push!(hosts, host)
		push!(pids, pid)
	end

    # remove worker processes
    for w in workers()
        rmprocs(w)
    end
end

main(ARGS)
