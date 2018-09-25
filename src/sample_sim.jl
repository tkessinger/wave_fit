#!/usr/bin/env julia

## sample_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Basic implementation of WaveFit.jl

using Revise
using Distributions, WaveFit

K = 1000 # carrying capacity

sigma = .05
delta = .01
s = .1
UL = 10
mu = 1e-4

landscape = Landscape(sigma, delta, s, UL, [mu,mu])
pop = Population(K,landscape)

tic()

for i in 1:100;
    evolve!(pop)
    if i%100 == 0
        println("$i $(length(pop.clones))")
    end
end

toc()

fcpop = FCPopulation(K,landscape)

Nhist = []
tic()

for i in 1:5000;
    evolve!(fcpop)
    push!(Nhist, sum([x.n for x in fcpop.classes]))
    if i%100 == 0
        println("$i $(length(fcpop.classes))")
    end
end

toc()
