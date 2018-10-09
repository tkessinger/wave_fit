#!/usr/bin/env julia

## WaveFit.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Simulation model for "traveling (fitness) wave" populations.


module WaveFit
export evolve!, selection!, mutation!
export Clone
export Population, Landscape
export FitnessClass, Population
export get_mean_fitness

using Distributions, Combinatorics

struct Landscape
    σ::Float64
    δ::Float64
    s::Float64
    β::Float64
    UL::Float64
    μ::Array{Float64,1}

    Landscape() = new(1e-6, 0.001, 0.01, 1.0, 10, [1e-4, 1e-4])
    Landscape(σ, δ, s, β, UL, μ) = new(σ, δ, s, β, UL, μ)
end

mutable struct FitnessClass
    n::Int64
    bg_mutations::Int64
    loci::Array{Int64,1}
    bg_fitness::Float64
    loci_fitness::Float64

    # constructors
    FitnessClass(n, bg_mutations, loci) = new(n, bg_mutations, loci, 0.0, 0.0)
    FitnessClass(n) = new(n, 0, [0, 0])
    FitnessClass() = FitnessClass(1)
end

function copy(fc::FitnessClass)
    return FitnessClass(fc.n, fc.bg_mutations, fc.loci)
end

# generate a key from a FitnessClass as an array
function key(fc::FitnessClass)
    return vcat([fc.bg_mutations], fc.loci)
end

mutable struct Population
    K::Int64    # carrying capacity
    N::Int64    # population size
    landscape::Landscape
    classes::Dict{Array{Int64,1}, FitnessClass}
    mean_bg_fitness::Float64
    mean_loci_fitness::Float64
    mean_fitness::Float64
    var_bg_fitness::Float64
    generation::Int64

    Population(K, landscape) = new(K, K, landscape, Dict(key(FitnessClass(K)) => FitnessClass(K)),
                                   1.0, 1.0, 1.0, 1.0, 0)
    Population(K) = Population(K, Landscape())
end

function loci_fitness(class::FitnessClass, pop::Population)
    # returns the sign-epistatic fitness for the focal loci
    return pop.landscape.s*(class.loci[1])*(class.loci[2]) -
            pop.landscape.δ*(class.loci[1])*(1 - class.loci[2]) -
            pop.landscape.δ*(1 - class.loci[1])*(class.loci[2])
end

function calc_fitness!(pop::Population)

    pop.mean_bg_fitness = 0
    pop.mean_loci_fitness = 0
    pop.mean_fitness = 0
    pop.var_bg_fitness = 0
    pop.N = sum([class.n for class in values(pop.classes)])

    # calculate the mean fitness and variance
    # there might be a smarter way to do this other than iterating over all clones twice...
    for class in values(pop.classes)
        class.bg_fitness = pop.landscape.β * class.bg_mutations
        class.loci_fitness = loci_fitness(class, pop)
        pop.mean_bg_fitness += class.bg_fitness * class.n / pop.N
        pop.var_bg_fitness += class.bg_fitness^2 * class.n / pop.N
        pop.mean_loci_fitness += class.loci_fitness * class.n / pop.N
    end
    pop.var_bg_fitness = pop.var_bg_fitness - pop.mean_bg_fitness^2

    if pop.var_bg_fitness > 0
        pop.mean_bg_fitness *= pop.landscape.σ / sqrt(pop.var_bg_fitness)
    end

    pop.mean_fitness = pop.mean_bg_fitness + pop.mean_loci_fitness
end

function selection!(pop::Population)
    # generates a Poisson-distributed offspring number for each class
    # still todo: incorporate the focal loci

    calc_fitness!(pop)

    # poisson offspring number distribution dependent on K and the fitness
    # the K/N term keeps the population size around K
    for (k, class) in pop.classes
        bg_fitness = class.bg_fitness
        if pop.var_bg_fitness > 0
           bg_fitness *= pop.landscape.σ/sqrt(var_bg_fitness)
        end
        total_fitness = bg_fitness + class.loci_fitness
        offspring_dist = Poisson(max(class.n * (1 + total_fitness - pop.mean_fitness) * pop.K / pop.N, 0.0))
        #offspring_dist = Poisson(max(class.n*(1+total_fitness)/(1+mean_fitness)*population.K/N,0.0))

        num_offspring = rand(offspring_dist)
        class.n = num_offspring

        # prune any empty fitness classes
        if class.n == 0
            pop!(pop.classes, k)
        end
    end
end

function mutation!(pop::Population)
    # randomly adds Poisson(UL) mutations to each individual.

    mut_dist = Poisson(pop.landscape.UL)
    # get the set of keys from the fitness class dict (`collect` is faster)
    ks = collect(keys(pop.classes))
    for k in ks
        class = pop.classes[k]
        tempclass = copy(class)
        for indv in 1:class.n;
            num_muts = rand(mut_dist)
            if num_muts > 0
                class.n -= 1
                tempclass.bg_mutations = class.bg_mutations + num_muts
                tempk = key(tempclass)
                if (tempk in keys(pop.classes))
                    pop.classes[tempk].n += 1
                else
                    pop.classes[tempk] = copy(tempclass)
                    pop.classes[tempk].n = 1
                end
            end
        end

        # class is empty now, so delete it
        if class.n == 0
            pop!(pop.classes, k)
        end
    end
end

function evolve!(pop::Population)
    selection!(pop)
    mutation!(pop)
    pop.generation += 1
end

# final end statement to close the module
end
