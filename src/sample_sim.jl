#!/usr/bin/env julia

## sample_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Basic implementation of WaveFit.jl

using Revise
using WaveFit
using Profile
using Distributions
using Stats

const K = 1000 # carrying capacity

const sigma = .05
const delta = .01
const s = .01
const UL = 100
const mu = 1e-4

landscape = Landscape(sigma, delta, s, UL, [mu,mu])
pop = Population(K,landscape)

@time (
for i in 1:10000;
    evolve!(pop)
end)

using Plots
bar([c.bg_mutations for c in values(pop.classes)], [c.n for c in values(pop.classes)])
