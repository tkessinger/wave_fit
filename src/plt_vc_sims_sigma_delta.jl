#!/usr/bin/env julia

## plt_vc_sims_weissman.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Reproduces main crossing time plots from Weissman et al. (2009).


using PyPlot, JLD, Glob, WaveFit, LaTeXStrings

results = Dict()

files = glob("output/vc_sims/dlx3/var_sigma_delta_*.jld2")
fixed_var_vals = []
fixed_vars = ["mu1", "mu2", "s", "K"]

indep_var = "delta"

for (fi, file) in enumerate(files)
    f = load(file)
    if fi == 1
        [push!(fixed_var_vals, f["parsed_args"][x]) for x in fixed_vars]
        #println("ok")
        #println(fixed_var_vals)
    end
    file_results = convert(Array{Int64,1},f["crossing_times"])
    indep_var_val = f["parsed_args"][indep_var]
    if indep_var_val == 0
        indep_var_val = 1e-8
    end
    tmp_sigma = f["parsed_args"]["sigma"]
    if length(f["crossing_times"]) < 50
        println(f["parsed_args"]["sigma"], "  ", f["parsed_args"]["delta"])
        println(length(f["crossing_times"]))
    end
    #println("$fig_to_make, $indep_var_val, $tmp_sigma, $(length(file_results))")
    #if f["params"]["sigma"] == sigmalist[1]
    #    push!(low_sigma_results, [indep_var_vals, file_results])
    # elseif f["params"]["sigma"] == sigmalist[2]
    #     push!(high_sigma_results, [indep_var_vals, file_results])
    if [tmp_sigma, indep_var_val] in collect(keys(results))
        [push!(results[[tmp_sigma, indep_var_val]], x) for x in file_results]
    else
        results[[tmp_sigma, indep_var_val]] = Float64[]
        [push!(results[[tmp_sigma, indep_var_val]], x) for x in file_results]

    end
end
#println(low_sigma_results)
#println(high_sigma_results)

indep_var_vals = unique([x[2] for x in collect(keys(results))])
sigma_vals = unique([x[1] for x in collect(keys(results))])
sort!(indep_var_vals)
sort!(sigma_vals)

fig = figure()

sigma_val_labels = [
L"1 \times 10^{-8}",
L"1 \times 10^{-2}",
L"5 \times 10^{-2}",
L"1 \times 10^{-1}",
L"5 \times 10^{-1}"
]


for (si, sigma) in enumerate(sigma_vals)
    #println([results[[sigma,x]] for x in indep_var_vals])
    asymmetric_err_length = [length(results[[sigma, x]]) for x in indep_var_vals]
    asymmetric_err_1 = floor.(Int64, 0.25*asymmetric_err_length)
    asymmetric_err_2 = floor.(Int64, 0.75*asymmetric_err_length)
    means = [median(results[[sigma, x]]) for x in indep_var_vals]
    sorted_results = [sort(results[[sigma, x]]) for x in indep_var_vals]
    #println(length(sorted_results))
    lower_err = [sorted_results[i][asymmetric_err_1[i]] for i in 1:length(sorted_results)]
    upper_err = [sorted_results[i][asymmetric_err_2[i]] for i in 1:length(sorted_results)]
    #println("made it here")
    #println([mean(results[[sigma, x]]) for x in indep_var_vals])
    #println([std(results[[sigma, x]]) for x in indep_var_vals])
    ax = errorbar(
    indep_var_vals,
    means,
    yerr = [means-lower_err,upper_err-means],
    label = sigma_val_labels[si]
    )
    println("$sigma, $means, $lower_err, $upper_err")
    println("$([means-lower_err, upper_err-means])")
end


#plot(indep_var_vals)

ylabel(L"\tau")
xlabel(L"\delta")
title("crossing time versus " * L"\delta")
xscale("log")
yscale("log")
legend(loc=2)
savefig("var_delta_sigma_v2.pdf")
#    return low_sigma_results, high_sigma_results

# figure()
# for (si, sigma) in enumerate(sigma_vals)
#     for (xi, x) in enumerate(indep_var_vals)
#         hist(results[[sigma,x]])
#     end
# end


# fig2 = figure()
# tmp_prediction = [pred_function(x, fixed_var_vals) for x in indep_var_vals]
# for sigma in sigma_vals
#     #sigma = minimum(sigma_vals)
#     means = [mean(results[[sigma, x]]) for x in indep_var_vals]
#     #plot(indep_var_vals, means)
#     ratio = [means[x]/tmp_prediction[x] for x in 1:(length(means))]
#
#     plot(indep_var_vals, ratio, label=sigma)
# end
# legend(loc=2)
# if (fig_to_make != 3)
#     xscale("log")
#     #yscale("log")
# end
# ylabel("ratio of simulation to analytical prediction")
# ylim([0,2])
# title("figure \($(collect('a':'z')[fig_to_make])\)")
# xlabel(indep_var_string)

#    return (sigmalist, UL_list, all_results)
