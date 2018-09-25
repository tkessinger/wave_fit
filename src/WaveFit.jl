#!/usr/bin/env julia

## WaveFit.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Simulation model for "traveling (fitness) wave" populations.

## STILL TO DO
## consider switching to "fitness-class" model

module WaveFit
export evolve!
export Clone
export Population, Landscape
export FitnessClass, FCPopulation

using Distributions, Combinatorics

struct Clone
    # new type for individual clones in population
    # each clone has a population size, genotype, and history
    # might be smart to make fitness an attribute as well, to avoid need to look up fitnesses constantly
    # struct is NOT mutable, be wary. we might change this
    # history is an array of 1x2 arrays: first element is the genotype, second is the mutation time

    n::Int64
    bg_mutations::Int64
    loci::Array{Int64,1}

    # constructors
    Clone() = new(1,0,[0,0])
    Clone(n) = new(n,0,[0,0])
    Clone(n,bg_mutations,loci) = new(n, bg_mutations, loci)
end

struct FitnessClass

    n::Int64
    bg_mutations::Int64
    loci::Array{Int64,1}

    # constructors
    FitnessClass() = new(1,0,[0,0])
    FitnessClass(n) = new(n,0,[0,0])
    FitnessClass(n,bg_mutations,loci) = new(n, bg_mutations, loci)

end

struct Landscape
    σ::Float64
    δ::Float64
    s::Float64
    UL::Float64
    μ::Array{Float64}

    Landscape() = new(1e-6,0.001,0.01,10,[1e-4,1e-4])
    Landscape(σ,δ,s,UL,μ) = new(σ,δ,s,UL,μ)
end

mutable struct FCPopulation
    #
    K::Int64
    landscape::Landscape
    classes::Array{FitnessClass,1}
    generation::Int64

    FCPopulation(K) = new(K,Landscape(),FitnessClass[FitnessClass(K)],0)
    FCPopulation(K,landscape) = new(K,landscape,FitnessClass[FitnessClass(K)],0)
end

mutable struct Population
    # type for storing both the population itself (clone list) and associated parameters
    K::Int64
    landscape::Landscape
    clones::Array{Clone,1}
    generation::Int64
    bubbles::Array{Int64,1}

    # constructor
    Population(K) = new(K,Landscape(),Clone[Clone(K)],0,Int64[0,0])
    Population(K,landscape) = new(K,landscape,Clone[Clone(K)],0,Int64[0,0])
end


function change_clone_size(clone::Clone, newsize::Int64)
    # copies an existing clone but modifies the size
    new_clone = Clone(newsize, clone.bg_mutations, clone.loci)
    return new_clone
end

function change_clone_size(class::FitnessClass, newsize::Int64)
    # copies an existing clone but modifies the size
    new_clone = FitnessClass(newsize, class.bg_mutations, class.loci)
    return new_clone
end


function loci_fitness(clone::Clone, population::Population)
    # returns the sign-epistatic fitness for the focal loci
    return population.landscape.δ*(clone.loci[1] != 0)*(clone.loci[2] == 0) + population.landscape.s*(clone.loci[1] != 0)*(clone.loci[2] != 0)
end


function selection!(population::FCPopulation)
    # generates a Poisson-distributed offspring number for each class
    # still todo: incorporate the focal loci

    new_classes = []
    mean_bg_fitness = 0
    mean_loci_fitness = 0
    mean_fitness = 0
    var_bg_fitness = 0
    N = sum([class.n for class in population.classes])

    # calculate the mean fitness and variance
    # there might be a smarter way to do this other than iterating over all clones twice...
    for (cli, class) in enumerate(population.classes);
        tmp_bg_fit = class.n*class.bg_mutations/N
        mean_bg_fitness += tmp_bg_fit
        var_bg_fitness += tmp_bg_fit^2/class.n
        #mean_loci_fitness += loci_fitness(class, population)*class.n/N
    end

    if var_bg_fitness > 0
        mean_bg_fitness *= population.landscape.σ^2/var_bg_fitness
    end

    mean_fitness = mean_bg_fitness + mean_loci_fitness

    # poisson offspring number distribution dependent on K and the fitness
    # the K/N term keeps the population size around K
    for (cli, class) in enumerate(population.classes);
        bg_fitness = 1.0*class.bg_mutations
        if var_bg_fitness > 0
           bg_fitness *= population.landscape.σ^2/var_bg_fitness
        end
        total_fitness = bg_fitness# + loci_fitness(clone, population)
        offspring_dist = Poisson(max(class.n*(1+total_fitness-mean_fitness)*population.K/N,0.0))
        num_offspring = rand(offspring_dist)

        tmp_clone = change_clone_size(class, num_offspring)
        # prune any empty clones
        if tmp_clone.n != 0
            push!(new_classes, tmp_clone)
        end
    end
    # update clone list
    population.classes = new_classes
end

function selection!(population::Population)
    # generates a Poisson-distributed offspring number for each clone
    # still todo: incorporate the focal loci

    new_clones = Clone[]
    mean_bg_fitness = 0
    mean_loci_fitness = 0
    mean_fitness = 0
    var_bg_fitness = 0
    N = sum([clone.n for clone in population.clones])

    # calculate the mean fitness and variance
    # there might be a smarter way to do this other than iterating over all clones twice...
    for (cli, clone) in enumerate(population.clones);
        tmp_bg_fit = clone.n*clone.bg_mutations/N
        mean_bg_fitness += tmp_bg_fit
        var_bg_fitness += tmp_bg_fit^2/clone.n
        mean_loci_fitness += loci_fitness(clone, population)*clone.n/N
    end

    if var_bg_fitness > 0
        mean_bg_fitness *= population.landscape.σ^2/var_bg_fitness
    end

    mean_fitness = mean_bg_fitness + mean_loci_fitness

    # poisson offspring number distribution dependent on K and the fitness
    # the K/N term keeps the population size around K
    for (cli, clone) in enumerate(population.clones);
        bg_fitness = 1.0*clone.bg_mutations
        if var_bg_fitness > 0
           bg_fitness *= population.landscape.σ^2/var_bg_fitness
        end
        total_fitness = bg_fitness + loci_fitness(clone, population)
        offspring_dist = Poisson(max(clone.n*(1+total_fitness-mean_fitness)*population.K/N,0.0))
        num_offspring = rand(offspring_dist)

        tmp_clone = change_clone_size(clone, num_offspring)
        # prune any empty clones
        if tmp_clone.n != 0
            push!(new_clones, tmp_clone)
        end
    end
    # update clone list
    population.clones = new_clones
end

function mutation!(population::Population)
    # randomly adds Poisson(UL) mutations to each individual.
    # there is a smarter way to do this.

    new_clones = Clone[]
    mut_dist = Poisson(population.landscape.UL)
    for (cli, clone) in enumerate(population.clones);
        parent_clone_size = clone.n
        for (ii, indv) in enumerate(1:clone.n);
            num_muts = rand(mut_dist)
            if num_muts > 0
                parent_clone_size -= 1
                push!(new_clones, Clone(1, clone.bg_mutations + num_muts, clone.loci))
            end
        end
        # delete this clone if it's empty
        if parent_clone_size != 0
            push!(new_clones, change_clone_size(clone, parent_clone_size))
        end
    end
    # update clone list
    population.clones = new_clones
end


function mutation!(population::FCPopulation)
    # randomly adds Poisson(UL) mutations to each individual.
    # there is a smarter way to do this.

    sort!(population.classes, by=x->x.bg_mutations)

    base_mutations = population.classes[1].bg_mutations
    new_classes = Int64[0]
    mut_dist = Poisson(population.landscape.UL)
    for (cli, class) in enumerate(reverse(population.classes));
        parent_clone_size = class.n
        for (ii, indv) in enumerate(1:class.n);
            num_muts = rand(mut_dist)
            if num_muts > 0
                parent_clone_size -= 1
                if class.bg_mutations + num_muts - base_mutations + 1 >=  length(new_classes)
                    #println("ok, $(class.bg_mutations), $num_muts, $base_mutations, $(length(new_classes))")
                    [push!(new_classes, 0) for x in 1:(length(new_classes)+class.bg_mutations+num_muts-base_mutations+1)]
                end
                #println("$(class.bg_mutations), $num_muts, $base_mutations, $(length(new_classes))")
                #println("$(class.bg_mutations+num_muts-base_mutations+1)")
                # if class.bg_mutations+num_muts-base_mutations+1 > length(new_classes)
                #     println("$(class.bg_mutations), $num_muts, $base_mutations, $(length(new_classes))")
                # end
                new_classes[class.bg_mutations+num_muts-base_mutations+1] += 1
                #push!(new_clones, Clone(1, clone.bg_mutations + num_muts, clone.loci))
            end
        end
        if class.bg_mutations + base_mutations + 1 >= length(new_classes)
            [push!(new_classes, 0) for x in 1:(length(new_classes)+class.bg_mutations-base_mutations+1)]
        end
        # if class.bg_mutations-base_mutations+1 > length(new_classes)
        #     println("$(class.bg_mutations), $base_mutations, $(length(new_classes))")
        # end
        new_classes[class.bg_mutations-base_mutations+1] += parent_clone_size
    end
    # update clone list
    population.classes = FitnessClass[]
    for (cli, class) in enumerate(new_classes);
        if class != 0
            push!(population.classes, FitnessClass(class, cli+base_mutations, [0,0]))
        end
    end
end


function cleanup!(population::Population)
    # idiotproofing function that removes empty clones
    for (cli, clone) in enumerate(population.clones);
        if clone.n == 0
            splice!(population.clones, cli)
        end
    end
end

function evolve!(population::Population)
    selection!(population)
    mutation!(population)
    population.generation += 1
end

function evolve!(population::FCPopulation)
    selection!(population)
    mutation!(population)
    population.generation += 1
end


# final end statement to close the module
end
