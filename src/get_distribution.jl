#!/usr/bin/env julia

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=30-00:00:00

# Set name of job shown in squeue
#SBATCH --job-name get_distribution.jl

# Request CPU resources
#SBATCH -p Long
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --account=8e300-jva238-lowprio

## get_distribution.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Gets and saves the size distribution and fixation probability of clones
## as a function of the background fitness.

#using Revise
using WaveFit
using Distributions
using ArgParse, JLD2
using Dates, Random

K = 1e4
sigma = 1e-1
delta = 0.0
s = 0.0
beta = 0.1
UL = 1.0

num_attempts = 10000

num_loci = 10

burn_time = 0.1*K

landscape = Landscape(sigma,
                      delta,
                      0.0,
                      beta,
                      UL,
                      [0.0, 0.0])

pop = Population(K, landscape)
for i in 1:burn_time
    evolve_multi!(pop)
end

# initialize mutations and run until valley crossing
segregating_clones = []
attempt_num = 1

clone_sizes = zeros(Int64, num_attempts)
clone_lifetimes = zeros(Int64, num_attempts)
clone_fitnesses = zeros(Float64, num_attempts)

clones_initiated = false

cleared_loci = Int64[]

fixed_loci = []

# we need to stop when the requisite number of attempts has been reached,
# but also continue tracking clones that already exist,
# even if that number has been reached
while attempt_num <= num_attempts || (length(segregating_clones) > 0 && clones_initiated)
    cutoff = min(num_loci, num_attempts - attempt_num + 1)
    # don't try to inject a new mutant if there are no wild type individuals
    while length(segregating_clones) < cutoff && sum([cl.n for cl in get_all_mutated_classes(pop)]) < pop.N
        global clones_initiated = true
        # mutate an individual and get the new clone's info
        tmp_clone_data = inject_mutation!(pop, attempt_num)
        new_clone_info = [attempt_num, tmp_clone_data[1], tmp_clone_data[2], pop.generation]
        push!(segregating_clones, new_clone_info)
        global attempt_num += 1
        println("$attempt_num")
    end
    # dummy array of clones to empty
    # purging an array while iterating over it would be a Very Bad Thing™
    clones_to_clear = []
    for (cl_id, cl_info) in enumerate(segregating_clones)
        k_id = cl_info[1]
        clone_size = sum([cl.n for cl in get_mutated_classes(pop, k_id)])
        clone_sizes[k_id] += clone_size
        # if a lineage has fixed, we need to discard it
        # and reset the population to wild type
        # store the fact that it fixed, though!
        if clone_size == pop.N
            push!(cleared_loci, k_id)
            fixed_clone_info = [k_id, cl_info[3], floor(Int64, pop.generation - cl_info[4])]
            push!(fixed_loci, fixed_clone_info)
            clear_loci!(pop, k_id)
            # reset all the clone's info
            tmp_clone_data = inject_mutation!(pop, k_id)
            new_clone_info = [k_id, tmp_clone_data[1], tmp_clone_data[2], pop.generation]
            segregating_clones[cl_id] = new_clone_info
            clone_sizes[k_id] = 1
        # if the clone is empty, we'll purge it later
        elseif clone_size == 0
            push!(clones_to_clear, cl_info)
            clone_lifetimes[k_id] = floor(Int64, pop.generation - cl_info[4])
            clone_fitnesses[k_id] = cl_info[3]
        end
    end
    # remove empty clones from segregating_clones
    for cl_info in clones_to_clear
        filter!(x->x≠cl_info,segregating_clones)
    end
    # increment the population
    evolve_multi!(pop)
end

for i in 1:num_attempts
   println("$(clone_sizes[i]), $(clone_lifetimes[i]), $(clone_fitnesses[i])")
end

empty_clones = []
for i in 1:num_attempts
   if clone_lifetimes[i] == 0
       push!(empty_clones, i)
   end
end
println("$empty_clones")
