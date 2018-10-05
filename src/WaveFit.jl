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
    UL::Float64
    μ::Array{Float64,1}

    Landscape() = new(1e-6,0.001,0.01,10,[1e-4,1e-4])
    Landscape(σ,δ,s,UL,μ) = new(σ,δ,s,UL,μ)
end

mutable struct FitnessClass

    n::Int64
    bg_mutations::Int64
    loci::Array{Int64,1}

    # constructors
    FitnessClass(n, bg_mutations, loci) = new(n, bg_mutations, loci)
    FitnessClass(n) = new(n,0,[0,0])
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
    #
    K::Int64
    landscape::Landscape
    classes::Dict{Array{Int64,1}, FitnessClass}
    generation::Int64

    Population(K,landscape) = new(K, landscape, Dict(key(FitnessClass(K)) => FitnessClass(K)), 0)
    Population(K) = Population(K, Landscape())
end

function get_mean_fitness(population::Population)

    mean_bg_fitness = 0
    mean_loci_fitness = 0
    mean_fitness = 0
    var_bg_fitness = 0
    N = sum([class.n for class in values(population.classes)])

    # calculate the mean fitness and variance
    # there might be a smarter way to do this other than iterating over all clones twice...
    for class in values(population.classes)
        tmp_bg_fit = class.n*class.bg_mutations/N
        mean_bg_fitness += tmp_bg_fit
        var_bg_fitness += tmp_bg_fit^2/class.n
        #mean_loci_fitness += loci_fitness(class, population)*class.n/N
    end
    if var_bg_fitness > 0
        mean_bg_fitness *= population.landscape.σ^2/var_bg_fitness
    end

    mean_fitness = mean_bg_fitness + mean_loci_fitness

    return [mean_bg_fitness, mean_fitness, var_bg_fitness]

end

function selection!(population::Population)
    # generates a Poisson-distributed offspring number for each class
    # still todo: incorporate the focal loci

    N = sum([class.n for class in values(population.classes)])
    (mean_bg_fitness, mean_fitness, var_bg_fitness) = get_mean_fitness(population)

    # poisson offspring number distribution dependent on K and the fitness
    # the K/N term keeps the population size around K
    for (k, class) in population.classes
        bg_fitness = 1.0*class.bg_mutations
        if var_bg_fitness > 0
           bg_fitness *= population.landscape.σ^2/var_bg_fitness
        end
        total_fitness = bg_fitness# + loci_fitness(clone, population)
        offspring_dist = Poisson(max(class.n*(1+total_fitness-mean_fitness)*population.K/N,0.0))
        #offspring_dist = Poisson(max(class.n*(1+total_fitness)/(1+mean_fitness)*population.K/N,0.0))

        num_offspring = rand(offspring_dist)
        class.n = num_offspring

        # prune any empty fitness classes
        if class.n == 0
            pop!(population.classes, k)
        end
    end
end

function mutation!(population::Population)
    # randomly adds Poisson(UL) mutations to each individual.
    # there is a smarter way to do this.

    mut_dist = Poisson(population.landscape.UL)
    # get the set of keys from the fitness class dict (NEED to `collect`)
    ks = collect(keys(population.classes))
    for k in ks
        class = population.classes[k]
        tempclass = copy(class)
        for indv in 1:class.n;
            num_muts = rand(mut_dist)
            if num_muts > 0
                class.n -= 1
                tempclass.bg_mutations = class.bg_mutations + num_muts
                tempk = key(tempclass)
                if (tempk in keys(population.classes))
                    population.classes[tempk].n += 1
                else
                    population.classes[tempk] = copy(tempclass)
                    population.classes[tempk].n = 1
                end
            end
        end

        # class is empty now, so delete it
        if class.n == 0
            pop!(population.classes, k)
        end
    end
end

function evolve!(population::Population)
    selection!(population)
    mutation!(population)
    population.generation += 1
end

# final end statement to close the module
end
