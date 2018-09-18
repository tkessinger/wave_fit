#!/usr/bin/env julia

## sample_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Basic implementation of WaveFit.jl

# tell Julia where the module is located
using Revise
using Distributions, WaveFit

K = 10000 # carrying capacity

sigma = .05
delta = .01
s = .1
UL = 10
mu = 1e-4

landscape = Landscape(sigma, delta, s, UL, [mu,mu])
pop = Population(K,landscape)

Nhist = []

tic()

for i in 1:10000;
    evolve!(pop)
    push!(Nhist, sum([x.n for x in pop.clones]))
    println("generation $i")
end

toc()
