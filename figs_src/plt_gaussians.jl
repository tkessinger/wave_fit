#!/usr/bin/env julia

## plt_vc_sims.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Open outputs of simulate_population.jl and plot results.
## This yields plots of crossing time versus Nσ for variable UL.

#ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"

using PyPlot, JLD, Glob, WaveFit, LaTeXStrings

function vec2array(v)
    r = length(v)
    c = length(v[1])
    #a = Array{Int64}(r,c)
    a = zeros(r,c)
    for i in 1:r, j in 1:c
        a[i,j] = v[i][j]
    end
    a
end

function get_mean_and_std(dist)
    #N = sum([dist[x,2] for x in 1:(size(dist)[1])])
    N = sum(dist[:,2])
    mean = 0
    var = 0
    [mean += dist[x,1]*dist[x,2]/N for x in 1:(size(dist)[1])]
    [var += dist[x,1]^2*dist[x,2]/N for x in 1:(size(dist)[1])]
    var -= mean^2
    return (mean, sqrt(var))
end

function normalize_and_sort(dist, sigma, mean, std)
    dist = convert.(Float64, dist)
    dist[:,1] -= mean
    dist[:,1] *= sigma/std
    dist[:,2] /= sum(dist[:,2])
    dist = sortrows(dist)
    return dist
end

date = "20190115"

    files = glob("output/vc_sims/" * date * "_sample_gaussian_*_*.jld2")
    # glob all the output files

    fitness_dists = Dict()
    fitness_means = Dict()
    fitness_devs = Dict()
    dist_lengths = Dict()
    fixed_dist_lengths = Dict()

    for file in files
        f = load(file)
        #println(f["params"])
        params = f["parsed_args"]
        sigma, UL, N, sampling_interval = params["sigma"], params["UL"], params["K"], params["sampling_interval"]
        if (sigma, UL, N) in keys(fitness_dists)
            transformed_dists = [vec2array(dist) for dist in f["pop_samples"]]
            [push!(fitness_dists[(sigma, UL, N)], x) for x in transformed_dists]
            means_and_vars = [get_mean_and_std(x) for x in transformed_dists]
            [push!(fitness_means[(sigma, UL, N)], x[1]) for x in means_and_vars]
            [push!(fitness_devs[(sigma, UL, N)], x[2]) for x in means_and_vars]
            #[push!(dist_lengths[(sigma, UL, N)], length(x)) for x in transformed_dists]
            push!(dist_lengths[(sigma, UL, N)], [length(x) for x in transformed_dists])

        else
            fitness_dists[(sigma, UL, N)] = []
            fitness_means[(sigma, UL, N)] = []
            fitness_devs[(sigma, UL, N)] = []
            dist_lengths[(sigma, UL, N)] = []
            transformed_dists = [vec2array(dist) for dist in f["pop_samples"]]
            [push!(fitness_dists[(sigma, UL, N)], x) for x in transformed_dists]
            means_and_vars = [get_mean_and_std(x) for x in transformed_dists]
            [push!(fitness_means[(sigma, UL, N)], x[1]) for x in means_and_vars]
            [push!(fitness_devs[(sigma, UL, N)], x[2]) for x in means_and_vars]
            #[push!(dist_lengths[(sigma, UL, N)], length(x)) for x in transformed_dists]
            push!(dist_lengths[(sigma, UL, N)], [length(x) for x in transformed_dists])
        end
        println("$file, $sigma, $UL, $N, $(length(f["pop_samples"]))")#", $(mean(f["pop_samples"]))")
    end
    sigmalist = sort(unique([x[1] for x in keys(fitness_dists)]))
    UL_list = sort(unique([x[2] for x in keys(fitness_dists)]))
    Nlist = sort(unique([x[3] for x in keys(fitness_dists)]))

    for (sigma, UL, N) in keys(fitness_dists)
        fixed_dist_lengths[(sigma, UL, N)] = []
        min_length = minimum([length(x) for x in dist_lengths[(sigma, UL, N)]])
        for i in 1:min_length
            tmp_lengths = dist_lengths[(sigma, UL, N)]
            push!(fixed_dist_lengths[(sigma, UL, N)], mean([tmp_lengths[1][i], tmp_lengths[2][i], tmp_lengths[3][i]]))
        end
    end

    #tau_grid = zeros(length(sigmalist), length(UL_list), length(Nlist))

    # for (si, sigma) in enumerate(sigmalist)
    #     for (ui, UL) in enumerate(UL_list)
    #         figure()
    #         for (Ni, N) in enumerate(Nlist)
    #             if (sigma, UL, N) in keys(fitness_dists)
    #                 title("$sigma, $UL")
    #                 tmp_dists = fitness_dists[(sigma,UL,N)]
    #                 tmp_means = fitness_means[(sigma,UL,N)]
    #                 tmp_devs = fitness_devs[(sigma,UL,N)]
    #                 for i in 1:length(tmp_dists)
    #                     if (i-1)%199 == 0
    #                         dist = normalize_and_sort(tmp_dists[i], sigma, tmp_means[i], tmp_devs[i])
    #                         plot(dist[:,1], dist[:,2], label = "N = "*string(N))
    #                         if any([x < 0 for x in dist[:,2]])
    #                             println("$sigma, $UL, $N, $i")
    #                         end
    #                     end
    #                 end
    #             end
    #             legend(loc=1)
    #         end
    #     end
    # end

    for (si, sigma) in enumerate(sigmalist)
        for (ui, UL) in enumerate(UL_list)
            if UL == 1.0
                figure()
                for (Ni, N) in enumerate(Nlist)
                    if (sigma, UL, N) in keys(fitness_dists)
                        title(latexstring("\$UL = $UL, \\sigma = 5 \\times 10^{$(round(Int64,log10(sigma/5)))}\$"))
                        plot(collect(1:length(fixed_dist_lengths[(sigma, UL, N)]))*.01, fixed_dist_lengths[(sigma, UL, N)], label=latexstring("\$N = 10^$(round(Int64,log10(N)))}\$"))
                    end
                end
                ylabel("number of classes")
                xlabel("generations/N")
                legend(loc=4)
                tight_layout()
                current_ylim = ylim()
                ylim([0, current_ylim[2]])
                vlines(.1, 0, current_ylim[2],linestyle = "--")
                savefig("equilibration_$(sigma)_$(UL).pdf")
            end
        end
    end


#     for (si, sigma) in enumerate(sigmalist)
#         for (ui, UL) in enumerate(UL_list)
#             figure()
#             for (Ni, N) in enumerate(Nlist)
#                 if (sigma, UL, N) in keys(fitness_dists)
#                     flattened_dist = zeros(2)'
#                     title("$sigma, $UL")
#                     tmp_dists = fitness_dists[(sigma,UL,N)]
#                     tmp_means = fitness_means[(sigma,UL,N)]
#                     tmp_devs = fitness_devs[(sigma,UL,N)]
#                     for i in 1:length(tmp_dists)
#                         #new_dist = normalize_and_sort(tmp_dists[i], sigma, tmp_means[i], tmp_devs[i])
# #                        for x in 1:(size(new_dist)[1])
# #                            flattened_dist = vcat(flattened_dist, new_dist[x,:]')
# #                        end
#                         flattened_dist = vcat(flattened_dist, tmp_dists[i])
#                     end
#                     #sortrows(flattened_dist)
#                 end
#                 #println("$flattened_dist")
#                 (mean, std) = get_mean_and_std(flattened_dist)
#                 flattened_dist = normalize_and_sort(flattened_dist, sigma, mean, std)
#                 plot(flattened_dist[:,1],flattened_dist[:,2])
#             end
#         end
#     end

    tmp_dists = fitness_dists[collect(keys(fitness_dists))[1]]
    tmp_means = fitness_means[collect(keys(fitness_dists))[1]]
    tmp_devs = fitness_devs[collect(keys(fitness_dists))[1]]
    i = 1
    #sigma = collect(keys(fitness_dists))[1][1]
    #dist = normalize_and_sort(tmp_dists[i], sigma, tmp_means[i], tmp_devs[i])


# figure()
# for (si, sigma) in enumerate(sigmalist)
#     for (ui, UL) in enumerate(UL_list)
#     ax = errorbar(Nlist, [mean(x) for x in tau_grid[si,ui,:]],
#         yerr = [std(x) for x in tau_grid[si,ui,:]],
#         label = latexstring("\$UL = $UL, \\sigma = 10^{$(round(Int64,log10(sigma)))}\$"))
#     end
# end
#
# ax = axes()
# legend(loc=2)
# ylabel(L"$\tau$")
# xlabel(L"$\sigma$")
# xscale("log")
# yscale("log")
# ylim([1e2,1e7])
# show()
# #savefig("neut_tunnel.pdf")
