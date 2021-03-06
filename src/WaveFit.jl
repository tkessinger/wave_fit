#!/usr/bin/env julia

## WaveFit.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Simulation model for "traveling (fitness) wave" populations.


module WaveFit
export evolve!, selection!, mutation!, evolve_multi!
export mutation_v2!, mutation_multinomial!, mutation_multinomial_v2!, mutation_multinomial_old!
export focal_mutation!
export Clone
export Population, Landscape
export FitnessClass
export get_mean_fitness, get_frequencies, get_mutated_classes, get_all_mutated_classes
export inject_mutation!, clear_loci!

using Distributions, StatsBase

mutable struct Landscape
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
    FitnessClass(n) = FitnessClass(n, 0, [0, 0])
    FitnessClass() = FitnessClass(1)
end

function Base.copy(fc::FitnessClass)
    return FitnessClass(fc.n, fc.bg_mutations, fc.loci)
end

function Base.copy!(fc_to::FitnessClass, fc_from::FitnessClass)
    fc_to.n = fc_from.n
    fc_to.bg_mutations = fc_from.bg_mutations
    fc_to.loci .= fc_from.loci
end

# generate a key from a FitnessClass as an array
# function key(fc::FitnessClass)
#     return vcat([fc.bg_mutations], fc.loci)
# end
#
# function key2(fc::FitnessClass)
#     return "$(fc.bg_mutations), $(fc.loci)"
# end
#
# function key3(fc::FitnessClass)
#     karr = zeros(Int64, 1+length(fc.loci))
#     karr[1] = fc.bg_mutations
#     karr[2:end] = fc.loci
#     return karr
# end

# borrowed from Steve Johnson
# https://discourse.julialang.org/t/poor-time-performance-on-dict/9656/12
struct FastHashUInt128; i::UInt128; end
Base.:(==)(x::FastHashUInt128, y::FastHashUInt128) = x.i == y.i
Base.hash(x::FastHashUInt128, h::UInt) = xor(x.i, h)

# generate a key from a FitnessClass as an Int128 assuming **biallelic** loci
function key(fc::FitnessClass)
    kint = UInt128(fc.bg_mutations)
    for i=1:length(fc.loci)
        @inbounds kint += (Int128(1) << (62+i)) * fc.loci[i]
    end
    return FastHashUInt128(kint)
end

# generate a key from an array with bg_mutations then loci
function key(arr::Array{Int64,1})
    kint = UInt128(arr[1])
    for i=1:length(arr)-1
        @inbounds kint += (Int128(1) << (62+i)) * arr[i+1]
    end
    return FastHashUInt128(kint)
end

mutable struct Population
    K::Int64    # carrying capacity
    N::Int64    # population size
    landscape::Landscape
    classes::Dict{FastHashUInt128, FitnessClass}
    mean_bg_fitness::Float64
    mean_loci_fitness::Float64
    mean_fitness::Float64
    var_bg_fitness::Float64
    generation::Int64

    Population(K, landscape) = new(K, K, landscape, Dict(key(FitnessClass(K)) => FitnessClass(K)),
                                   1.0, 1.0, 1.0, 1.0, 1)
    Population(K) = Population(K, Landscape())
end

function calc_pop_size!(pop::Population)
    pop.N = sum([class.n for class in values(pop.classes)])
end

function loci_fitness(class::FitnessClass, pop::Population)
    # returns the sign-epistatic fitness for the focal loci
    return pop.landscape.s*(1-(class.loci[1] != 0))*(1-(class.loci[2] != 0)) -
            pop.landscape.δ*(1 -(class.loci[1] != 0))*(class.loci[2] != 0) -
            pop.landscape.δ*(class.loci[1] != 0)*(1 - (class.loci[2] != 0))
end

function calc_fitness!(pop::Population)

    pop.mean_bg_fitness = 0
    pop.mean_loci_fitness = 0
    pop.mean_fitness = 0
    pop.var_bg_fitness = 0

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

    calc_pop_size!(pop)
    calc_fitness!(pop)
    # poisson offspring number distribution dependent on K and the fitness
    # the K/N term keeps the population size around K
    for (k, class) in pop.classes
        bg_fitness = class.bg_fitness
        if pop.var_bg_fitness > 0
           bg_fitness *= pop.landscape.σ/sqrt(pop.var_bg_fitness)
        end
        total_fitness = bg_fitness + class.loci_fitness
        offspring_dist = Poisson(max(class.n * (1 + total_fitness - pop.mean_fitness) * pop.K / pop.N, 0.0))

        num_offspring = rand(offspring_dist)
        class.n = num_offspring

        # prune any empty fitness classes
        if class.n == 0
            pop!(pop.classes, k)
        end
    end
    # population size needs to be updated after selection step
    calc_pop_size!(pop)
end

function mutation!(pop::Population)
    # randomly adds Poisson(UL) mutations to each individual.

    mut_dist = Poisson(pop.landscape.UL)

    # capture class keys and their sizes and iterate over each class
    keys_ns = [(k, pop.classes[k].n) for k in keys(pop.classes)]
    for (k, n) in keys_ns
        class = pop.classes[k]
        tempclass = copy(class)
        for indv in 1:n
            num_muts = rand(mut_dist)
            if num_muts > 0
                class.n -= 1
                tempclass.bg_mutations = class.bg_mutations + num_muts
                tempk = key(tempclass)
                if (tempk in keys(pop.classes))
                    pop.classes[tempk].n += 1
                else
                    # the `copy` here is important
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

function mutation_v2!(pop::Population)
    # randomly adds Poisson(UL) mutations to each individual.

    mut_dist = Poisson(pop.landscape.UL)
    num_muts = rand(mut_dist, pop.N)

    # capture class keys and their sizes and iterate over each class
    keys_ns = [(k, pop.classes[k].n) for k in keys(pop.classes)]
    i = 0
    for (k, n) in keys_ns
        class = pop.classes[k]
        tempclass = copy(class)
        for indv in 1:n
            i += 1
            if num_muts[i] > 0
                class.n -= 1
                tempclass.bg_mutations = class.bg_mutations + num_muts[i]
                tempk = key(tempclass)
                if (tempk in keys(pop.classes))
                    pop.classes[tempk].n += 1
                else
                    # the `copy` here is important
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

function mutation_multinomial!(pop::Population)
    # randomly adds Poisson(UL) mutations to each individual

    mut_dist = Poisson(pop.N*pop.landscape.UL)
    total_muts = rand(mut_dist)

    # create multinomial sample
    # (empirical evidence suggests this might be faster than
    #  Multinomial(n,k) for 10n<=k)
    distN = DiscreteUniform(1, pop.N)
    num_muts = zeros(Int64, pop.N)
    for i = 1:total_muts
        @inbounds num_muts[rand(distN)] += 1
    end

    # capture class keys and their sizes and iterate over each class
    keys_ns = [(k, pop.classes[k].n) for k in keys(pop.classes)]

    # iterate over each individual and update its fitness class
    i = 0
    for (k, n) in keys_ns
        class = pop.classes[k]
        tempclass = copy(class)
        for indv in 1:n
            i += 1
            if num_muts[i] > 0
                class.n -= 1
                tempclass.bg_mutations = class.bg_mutations + num_muts[i]
                tempk = key(tempclass)
                oldclass = get(pop.classes, tempk, nothing)
                if oldclass != nothing
                    oldclass.n += 1
                else # create new class with `copy`
                    newclass = copy(tempclass)
                    newclass.n = 1
                    pop.classes[tempk] = newclass
                end
            end
        end

        # class is empty now, so delete it
        if class.n == 0
            pop!(pop.classes, k)
        end
    end
end

function mutation_multinomial_v2!(pop::Population)
    # randomly adds Poisson(UL) mutations to each individual

    mut_dist = Poisson(pop.N*pop.landscape.UL)
    total_muts = rand(mut_dist)

    # create multinomial sample
    num_muts = rand(Multinomial(total_muts, pop.N))

    # capture class keys and their sizes and iterate over each class
    keys_ns = [(k, pop.classes[k].n) for k in keys(pop.classes)]

    # iterate over each individual and update its fitness class
    i = 0
    for (k, n) in keys_ns
        class = pop.classes[k]
        tempclass = copy(class)
        for indv in 1:n
            i += 1
            if num_muts[i] > 0
                class.n -= 1
                tempclass.bg_mutations = class.bg_mutations + num_muts[i]
                tempk = key(tempclass)
                if (tempk in keys(pop.classes))
                    pop.classes[tempk].n += 1
                else
                    # the `copy` here is important
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

function mutation_multinomial_old!(pop::Population, fixed_first_locus::Bool=false)
    # randomly adds Poisson(UL) mutations to each individual.

    mut_dist = Poisson(pop.N*pop.landscape.UL)
    num_muts = rand(mut_dist)

    muts_to_assign_dist = Multinomial(num_muts, [x.n for x in values(pop.classes)]/pop.N)
    muts_to_assign = rand(muts_to_assign_dist)

    # capture class keys and their sizes and iterate over each class
    kis_keys_ns = [(ki, k, pop.classes[k].n) for (ki, k) in enumerate(keys(pop.classes))]
    for (ki, k, n) in kis_keys_ns
        class = pop.classes[k]
        tempclass = copy(class)
        tmp_muts = muts_to_assign[ki]
        tmp_mut_dist = Multinomial(tmp_muts, n)
        class_muts = rand(tmp_mut_dist)
        for i in 1:n
            if class_muts[i] > 0
                class.n -= 1
                tempclass.bg_mutations = class.bg_mutations + class_muts[i]
                tempk = key(tempclass)
                if (tempk in keys(pop.classes))
                    pop.classes[tempk].n += 1
                else
                    # the `copy` here is important
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

function get_frequencies(pop::Population)
    freqs = Float64[0.0,0.0]
    for (k, class) in pop.classes
        freqs += [1.0*(x != false) for x in class.loci]*class.n
    end
    freqs /= pop.N
    return freqs
end

function get_mutated_classes(pop::Population, k_id::Int64)
    mutated_classes = FitnessClass[]
    locus = 1
    for (k, class) in pop.classes
        if class.loci[locus] == k_id
            push!(mutated_classes, class)
        end
    end
    return mutated_classes
end

function get_all_mutated_classes(pop::Population)
    mutated_classes = FitnessClass[]
    locus = 1
    for (k, class) in pop.classes
        if class.loci[locus] != 0
            push!(mutated_classes, class)
        end
    end
    return mutated_classes
end

function inject_mutation!(pop::Population, k_id::Int64=1)
    wt_muts = 1
    wt_classes = FitnessClass[]
    for (k, class) in pop.classes
        if class.loci[1] == 0
            push!(wt_classes, class)
        end
    end
    class_bg_fitness = 0.0
    mutated_class = []
    if length(wt_classes) > 0
        probabilities = weights(1.0*[x.n for x in wt_classes])
        classes_to_mutate = sample(wt_classes, probabilities, wt_muts)
        for cl in classes_to_mutate
            k = key([cl.bg_mutations, k_id, 0])
            oldclass = get(pop.classes, k, nothing)
            if oldclass != nothing
                oldclass.n += 1
            else
                tempclass = copy(cl)
                tempclass.n = 1
                tempclass.loci = [k_id, 0]
                pop.classes[k] = tempclass
                mutated_class = tempclass
            end
            calc_fitness!(pop)
            bg_fitness = pop.classes[k].bg_fitness
            if pop.var_bg_fitness > 0
               bg_fitness *= pop.landscape.σ/sqrt(pop.var_bg_fitness)
            end
            class_bg_fitness = bg_fitness - pop.mean_fitness
            cl.n -= 1
        end
    end
    return mutated_class, class_bg_fitness
end

function clear_loci!(pop::Population, k_id::Int64=1)
    classes_to_clear = []
    for (k, class) in pop.classes
        if class.loci[1] == k_id
            push!(classes_to_clear, class)
        end
    end
    for cl in classes_to_clear
        k = key([cl.bg_mutations, 0, 0])
        oldclass = get(pop.classes, k, nothing)
        if oldclass != nothing
            oldclass.n += cl.n
        else
            tempclass = copy(cl)
            tempclass.n = cl.n
            tempclass.loci = [0, 0]
            pop.classes[k] = tempclass
        end
        cl.n = 0
    end
end

function clear_loci!(pop::Population, cl::FitnessClass)

    k = key([cl.bg_mutations, 0, 0])
    oldclass = get(pop.classes, k, nothing)
    if oldclass != nothing
        oldclass.n += cl.n
    else
        tempclass = copy(cl)
        tempclass.n = cl.n
        tempclass.loci = [0, 0]
        pop.classes[k] = tempclass
    end
    cl.n = 0
end

function focal_mutation!(pop::Population, fixed_first_locus::Bool=false)
    current_freqs = get_frequencies(pop)
    wild_type = (1.0-current_freqs[1])*pop.N
    if wild_type < 0
        println(current_freqs)
        println(wild_type)
    end

    if !(fixed_first_locus)
        wt_dist = Poisson(wild_type*pop.landscape.μ[1])
        wt_muts = rand(wt_dist)
    elseif fixed_first_locus && get_frequencies(pop)[1] == 0
        wt_muts = 1
    elseif fixed_first_locus && get_frequencies(pop)[1] > 0
        wt_muts = 0
    end

    wt_classes = FitnessClass[]
    for (k, class) in pop.classes
        if class.loci[1] == 0
            push!(wt_classes, class)
        end
    end
    if length(wt_classes) > 0
        probabilities = weights(1.0*[x.n for x in wt_classes])
        classes_to_mutate = sample(wt_classes, probabilities, wt_muts)
        for cl in classes_to_mutate
            k = key([cl.bg_mutations, 1, 0])
            oldclass = get(pop.classes, k, nothing)
            if oldclass != nothing
                oldclass.n += 1
            else
                tempclass = copy(cl)
                tempclass.n = 1
                tempclass.loci = [1, 0]
                pop.classes[k] = tempclass
            end
            cl.n -= 1
        end
    end

    current_freqs = get_frequencies(pop)
    single_mut = (current_freqs[1]-current_freqs[2])*pop.N
    if single_mut < 0
        println(current_freqs)
        println(single_mut)
    end

    sm_dist = Poisson(single_mut*pop.landscape.μ[2])
    sm_muts = rand(sm_dist)
    sm_classes = FitnessClass[]
    for (k, class) in pop.classes
        if class.loci[2] == 0 && class.loci[1] != 0
            push!(sm_classes, class)
        end
    end
    if length(sm_classes) > 0
        probabilities = weights(1.0*[x.n for x in sm_classes])
        classes_to_mutate = sample(sm_classes, probabilities, sm_muts)
        for cl in classes_to_mutate
            k = key([cl.bg_mutations, 1, 1])
            oldclass = get(pop.classes, k, nothing)
            if oldclass != nothing
                oldclass.n += 1
            else
                tempclass = copy(cl)
                tempclass.n = 1
                tempclass.loci = [1, 1]
                pop.classes[k] = tempclass
            end
            cl.n -= 1
        end
    end
end

function evolve!(pop::Population)
    selection!(pop)
    mutation!(pop)
    focal_mutation!(pop)
    pop.generation += 1
end

function evolve_multi!(pop::Population, fixed_first_locus::Bool=false)
    selection!(pop)
    mutation_multinomial!(pop)
    focal_mutation!(pop, fixed_first_locus)
    pop.generation += 1
end

# final end statement to close the module
end
