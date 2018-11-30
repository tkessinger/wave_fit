#!/usr/bin/env julia

## sample_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Basic implementation of WaveFit.jl

using Revise
using WaveFit
#using Profile
using Distributions
using Stats

const K = 10000 # carrying capacity

const sigma = .05
const delta = .001
const s = .01
const UL = 10
const beta = 0.01
const mu = 0.0

numgens = 100000

landscape = Landscape(sigma, delta, s, UL, beta, [mu,mu])
pop = Population(K,landscape)

@time (
#while get_frequencies(pop)[2] < 0.5
for i in 1:numgens;
    evolve!(pop)
    if (pop.generation%(numgens/10) == 0)
        #bar([c.bg_mutations for c in values(pop.classes)], [c.n for c in values(pop.classes)],yscale=:log10)
        #bg_muts = [c.bg_mutations for c in values(pop.classes)]
        #println(maximum(bg_muts) - minimum(bg_muts))
        println(pop.generation)
    end
end
)
println("$K $numgens")

landscape = Landscape(sigma, delta, s, UL, beta, [mu,mu])
pop = Population(K,landscape)

@time (
#while get_frequencies(pop)[2] < 0.5
for i in 1:numgens;
    evolve_multi!(pop)
    if (pop.generation%(numgens/10) == 0)
        #bar([c.bg_mutations for c in values(pop.classes)], [c.n for c in values(pop.classes)],yscale=:log10)
        #bg_muts = [c.bg_mutations for c in values(pop.classes)]
        #println(maximum(bg_muts) - minimum(bg_muts))
        println(pop.generation)
    end
end
)
println("$K $numgens")


#using Plots
# plot the distribution of background mutations
#bar([c.bg_mutations for c in values(pop.classes)], [c.n for c in values(pop.classes)],yscale=:log10)
