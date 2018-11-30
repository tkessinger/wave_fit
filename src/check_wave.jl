#!/usr/bin/env julia

## simulate_population.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Check to make sure the mean fitness is advancing at a rate σ^2.

#using Revise
using WaveFit
using Distributions
using Stats
using ArgParse, JLD

using PyPlot

K = 1000000
sigma, delta, s, UL, beta = 1e-2, 1e-2, 1e-2, 100, 1e-2

mu1, mu2 = 0.0, 0.0

#burn_time = K/10
burn_time = 0

landscape = Landscape(sigma, delta, s, UL, beta, [0.0, 0.0])
pop = Population(K,landscape)

crossing_times = []

mean_fitness = []

while pop.generation < 10000
    evolve!(pop)
    if pop.generation > burn_time
        push!(mean_fitness, [pop.mean_bg_fitness, pop.mean_fitness])
    end
end

plot([x[1] for x in mean_fitness])
plot([x[2] for x in mean_fitness])
