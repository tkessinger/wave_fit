using StatsKit
using Revise
using BenchmarkTools
using Plots
using WaveFit

function test_randn(N::Int64, n::Int64)
    x = zeros(N)
    for i = 1:n
        x[rand(1:N)] += 1
    end

    return x
end

function test_alias(N::Int64, n::Int64)
    a = Distributions.AliasTable(fill(1/N, N))
    x = zeros(N)
    for i = 1:n
        x[rand(a)] += 1
    end
end

function test_mult_bare(N::Int64, n::Int64)
    x = zeros(N)
    Distributions.multinom_rand!(n, fill(1/N, N), x)

    return x
end

function test_mult_user(N::Int64, n::Int64)
    return rand(Multinomial(n, N))
end

# @btime test_randn(10000, 10000)
# @btime test_alias(10000, 10000)
# @btime test_mult_bare(10000, 10000)
# @btime test_mult_user(10000, 10000)

function test_r_poisson(λ::Real, n::Int64)
    return rand(Poisson(λ), n)
end

# maybe less accuracy than rpois?
function test_alias_poisson(λ::Real, n::Int64)
    probs = pdf.(Poisson(λ), 1:quantile(Poisson(λ),1-2*eps(1.0)))
    probs /= sum(probs)
    s = sampler(Categorical(probs))

    return rand(s, n)
end

# @btime test_r_poisson(100, 1000000)
# @btime test_alias_poisson(100, 1000000)

function grid_bench(Nvals, nvals, func)
    times = zeros(d,d)
    for i = 1:length(Nvals)
        for j = 1:length(nvals)
            println("benching $func at $(Nvals[i]), $(nvals[j])")
            bm = @benchmarkable $func($Nvals[$i], $nvals[$j])
            t = run(bm, samples=100)
            times[i,j] = median(t).time
        end
    end

    return times
end

d = 10
Nvals = nvals = round.(Int64, 10 .^ range(1, stop=6, length=d))

randn_times = grid_bench(Nvals, nvals, test_randn)
alias_times = grid_bench(Nvals, nvals, test_alias)
mult_bare_times = grid_bench(Nvals, nvals, test_mult_bare)
mult_user_times = grid_bench(Nvals, nvals, test_mult_user)

# ratio of the two multinomial tests. they look very similar though mult_bare slightly better
heatmap(
    mult_user_times./mult_bare_times,
    aspect_ratio=1,
    ticks=(range(1, stop=10, length=6), [10, 100, 1000,10000, 100000, 1000000]),
    xlabel="n", ylabel="k",
    clims=(0.9, 1.2),
    c=cgrad(:RdBu_r, :colorbrewer,
        unique(vcat(range(0.0, stop=0.333, length=6), range(0.333, stop=1.0, length=6))))
    )

# the alias method is always slower (not surprising), particularly for many samples
heatmap(
    log10.(alias_times./randn_times),
    aspect_ratio=1,
    ticks=(range(1, stop=10, length=6), [10, 100, 1000,10000, 100000, 1000000]),
    xlabel="n", ylabel="k",
    clims=(0,1.2),
    c=cgrad(:Reds, :colorbrewer)
    )

# randn multinomial tends to be faster than multinomial when 10n<=k
heatmap(
    log10.(randn_times./mult_bare_times),
    aspect_ratio=1,
    ticks=(range(1, stop=10, length=6), [10, 100, 1000,10000, 100000, 1000000]),
    xlabel="n", ylabel="k",
    c=cgrad(:RdBu_r, :colorbrewer,
        unique(vcat(range(0.0, stop=0.333, length=6), range(0.333, stop=1.0, length=6)))),
    clims=(-2,4)
    )

contour(
    log10.(randn_times./mult_bare_times),
    aspect_ratio=1,
    ticks=(range(1, stop=10, length=6), [10, 100, 1000,10000, 100000, 1000000]),
    xlabel="n", ylabel="k",
    c=cgrad(:RdBu_r, :colorbrewer,
        unique(vcat(range(0.0, stop=0.333, length=6), range(0.333, stop=1.0, length=6)))),
    clims=(-2,4), fill=true)

heatmap(
    log10.(randn_times),
    aspect_ratio=1,
    ticks=(range(1, stop=10, length=6), [10, 100, 1000,10000, 100000, 1000000]),
    xlabel="n", ylabel="k"
    )

heatmap(
    log10.(mult_bare_times),
    aspect_ratio=1,
    ticks=(range(1, stop=10, length=6), [10, 100, 1000,10000, 100000, 1000000]),
    xlabel="n", ylabel="k"
    )


function run_mutation(N, landscape, iters, mfunc)
    pop = Population(N, landscape)
    for i in 1:iters
        mfunc(pop)
    end
end


N = 1000
iters = 1000

# UL = 0.1
landscape = Landscape(1.0, 0.001, 0.1, 0.001, 0.1, [0.0, 0.0])
println("UL = 0.1")
print("mutation!                 ")
@btime run_mutation(N, landscape, iters, mutation!)
print("mutation_v2!              ")
@btime run_mutation(N, landscape, iters, mutation_v2!)
print("mutation_multinomial!     ")
@btime run_mutation(N, landscape, iters, mutation_multinomial!)
print("mutation_multinomial_v2!  ")
@btime run_mutation(N, landscape, iters, mutation_multinomial_v2!)
print("mutation_multinomial_old! ")
@btime run_mutation(N, landscape, iters, mutation_multinomial_old!)

# UL = 1
landscape = Landscape(1.0, 0.001, 0.1, 0.001, 1, [0.0, 0.0])
println("UL = 1")
print("mutation!                 ")
@btime run_mutation(N, landscape, iters, mutation!)
print("mutation_v2!              ")
@btime run_mutation(N, landscape, iters, mutation_v2!)
print("mutation_multinomial!     ")
@btime run_mutation(N, landscape, iters, mutation_multinomial!)
print("mutation_multinomial_v2!  ")
@btime run_mutation(N, landscape, iters, mutation_multinomial_v2!)
print("mutation_multinomial_old! ")
@btime run_mutation(N, landscape, iters, mutation_multinomial_old!)

# UL = 10
landscape = Landscape(1.0, 0.001, 0.1, 0.001, 10, [0.0, 0.0])
println("UL = 10")
print("mutation!                 ")
@btime run_mutation(N, landscape, iters, mutation!)
print("mutation_v2!              ")
@btime run_mutation(N, landscape, iters, mutation_v2!)
print("mutation_multinomial!     ")
@btime run_mutation(N, landscape, iters, mutation_multinomial!)
print("mutation_multinomial_v2!  ")
@btime run_mutation(N, landscape, iters, mutation_multinomial_v2!)
print("mutation_multinomial_old! ")
@btime run_mutation(N, landscape, iters, mutation_multinomial_old!)
