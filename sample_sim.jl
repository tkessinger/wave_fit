#!/usr/bin/env julia

## sample_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Basic implementation of WaveFit.jl

using Revise
#using Profile, ProfileView
using ProfileView
using Distributions, WaveFit

const K = 1000 # carrying capacity

const sigma = .05
const delta = .01
const s = .1
const UL = 100
const mu = 1e-4

landscape = Landscape(sigma, delta, s, UL, [mu,mu])
pop = Population(K,landscape)

# tic()
#
# for i in 1:10;
#     evolve!(pop)
#     if i%100 == 0
#         println("$i $(length(pop.clones))")
#     end
# end
#
# toc()

fcpop = FCPopulation(K,landscape)

Nhist = []
#tic()
println("list method")
@profile (for j = 1:2;
for i in 1:1000;
    evolve_list!(fcpop)
    #push!(Nhist, sum([x.n for x in fcpop.classes]))
    #if i%100 == 0
    #    println("$i $(length(fcpop.classes)) $(sum([x.n for x in fcpop.classes])), $(get_mean_fitness(fcpop))")
    #end
end
end)
println("list method profile")
Profile.print()

#toc()

fcpop = FCPopulation(K,landscape)

Nhist = []
#tic()

println("dict method")
@profile (for j = 1:2;
for i in 1:1000;
    evolve_dict!(fcpop)
    #push!(Nhist, sum([x.n for x in fcpop.classes]))
    #if i%100 == 0
    #    println("$i $(length(fcpop.classes)) $(sum([x.n for x in fcpop.classes])), $(get_mean_fitness(fcpop))")
    #end
end
end)
println("dict method profile")
Profile.print()
#toc()
