#!/usr/bin/env julia

## plt_vc_sims_weissman.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Reproduces main crossing time plots from Weissman et al. (2009).


using PyPlot, JLD, Glob, WaveFit, LaTeXStrings

function weissman_1(N, fixed_params)
    mu1, mu2, s, delta = fixed_params[1], fixed_params[2], fixed_params[3], fixed_params[4]
    println("N = $N, mu1 = $mu1, mu2 = $mu2, s = $s, delta = $delta")
    if N < 1/sqrt(mu1*s)
        return 1.0/mu1 + 1.0/(N*mu1*(s+delta))
    elseif N < 1/mu1
        return 1.0/(N*mu1*sqrt(mu2*s))
    else
        return sqrt(pi/(2.0*N*mu1*mu2*s))
    end
end

function weissman_2(N, fixed_params)
    mu1, mu2, s, delta = fixed_params[1], fixed_params[2], fixed_params[3], fixed_params[4]
    println("N = $N, mu1 = $mu1, mu2 = $mu2, s = $s, delta = $delta")
    if N < 1/delta*log(1+(delta*(exp(delta)-1)/(mu2*s)))
        return 1.0/(N*mu1*(exp(delta)-1)/(exp(N*delta)-1)) + 1.0/(N*mu1*(s+delta))
        #return 1.0/(N*mu1*(exp(delta)-1)/(exp(N*delta)-1))
    elseif N < 2*delta^2/(pi*mu1*mu2*s)
        return delta/(N*mu1*mu2*s)
    else
        return sqrt(pi/(2.0*N*mu1*mu2*s))
    end
end

function weissman_3(delta, fixed_params)
    N, mu1, mu2, s = fixed_params[1], fixed_params[2], fixed_params[3], fixed_params[4]
    println("N = $N, mu1 = $mu1, mu2 = $mu2, s = $s, delta = $delta")
    p = (-delta + sqrt(delta^2 + 4*mu2*s))/2
    return 1.0/(N*mu1*p)
    # if delta < 2*sqrt(mu2*s)
    #     return 1.0/(N*mu1*sqrt(mu2*s))
    # else
    #     return delta/(N*mu1*mu2*s)
    # end
end

function weissman_4(mu2, fixed_params)
    N, mu1, s, delta = fixed_params[1], fixed_params[2], fixed_params[3], fixed_params[4]
    println("N = $N, mu1 = $mu1, mu2 = $mu2, s = $s, delta = $delta")
#    return delta/(N*mu1*mu2*s)
    p = (-delta + sqrt(delta^2 + 4*mu2*s))/2
    return 1.0/(N*mu1*p)


end

function load_and_plot(fig_to_make::Int64)
    if fig_to_make==1
        indep_var = "K"
        indep_var_string = L"N"
        fixed_vars = ["mu1","mu2","s","delta"]
        pred_function = weissman_1
    elseif fig_to_make==2
        indep_var = "K"
        indep_var_string = L"N"
        fixed_vars = ["mu1","mu2","s","delta"]
        pred_function = weissman_2
    elseif fig_to_make==3
        indep_var = "delta"
        indep_var_string = L"\delta"
        fixed_vars = ["K","mu1","mu2","s"]
        pred_function = weissman_3
    elseif fig_to_make==4
        indep_var = "mu2"
        indep_var_string = L"\mu_2"
        fixed_vars = ["K","mu1","s","delta"]
        pred_function = weissman_4
    end
    #sigmalist = [1e-8,1e-2]
    #low_sigma_results = []
    #high_sigma_results = []

    results = Dict()

    namestring = "weissman_fig_"*string(fig_to_make)
    files = glob("output/vc_sims/" * namestring * "*.jld2")
    # glob all the output files

    fixed_var_vals = []

    for (fi, file) in enumerate(files)
        f = load(file)
        if fi == 1
            [push!(fixed_var_vals, f["parsed_args"][x]) for x in fixed_vars]
            println("ok")
            println(fixed_var_vals)
        end
        file_results = convert(Array{Int64,1},f["crossing_times"])
        indep_var_val = f["parsed_args"][indep_var]
        tmp_sigma = f["parsed_args"]["sigma"]
        println("$fig_to_make, $indep_var_val, $tmp_sigma, $(length(file_results))")
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

    println("made it here")

    # fig = figure()
    # ax1 = errorbar(
    #     [x[1] for x in low_sigma_results],
    #     [mean(x[2]) for x in low_sigma_results],
    #     yerr = [std(x[2]) for x in low_sigma_results],
    #     label = L"\sigma = 10^{-8}"
    # )
    # ax2 = errorbar(
    #     [x[1] for x in high_sigma_results],
    #     [mean(x[2]) for x in high_sigma_results],
    #     yerr = [std(x[2]) for x in high_sigma_results],
    #     label = L"\sigma = 10^{-2}"
    # )

    fig = figure()
    for sigma in sigma_vals
        #println([results[[sigma,x]] for x in indep_var_vals])
        asymmetric_err_length = [length(results[[sigma, x]]) for x in indep_var_vals]
        asymmetric_err_1 = floor.(Int64, 0.25*asymmetric_err_length)
        asymmetric_err_2 = floor.(Int64, 0.75*asymmetric_err_length)
        means = [mean(results[[sigma, x]]) for x in indep_var_vals]
        sorted_results = [sort(results[[sigma, x]]) for x in indep_var_vals]
        println(length(sorted_results))
        lower_err = [sorted_results[i][asymmetric_err_1[i]] for i in 1:length(sorted_results)]
        upper_err = [sorted_results[i][asymmetric_err_2[i]] for i in 1:length(sorted_results)]
        #println("made it here")
        #println([mean(results[[sigma, x]]) for x in indep_var_vals])
        #println([std(results[[sigma, x]]) for x in indep_var_vals])
        ax = errorbar(
        indep_var_vals,
        means,
        yerr = [means-lower_err,upper_err-means],
        label = sigma
        )
    end


    if fig_to_make==3
        plot_line_vals = collect(linspace(minimum(indep_var_vals),maximum(indep_var_vals),100))
    else
        plot_line_vals = collect(logspace(log10(minimum(indep_var_vals)),
            log10(maximum(indep_var_vals)),100))
    end
    println("$plot_line_vals")


    analytical_prediction = [pred_function(x, fixed_var_vals) for x in plot_line_vals]
    plot(plot_line_vals,analytical_prediction,c="k",ls="--")

    #plot(indep_var_vals)

    ylabel(L"\tau")
    xlabel(indep_var_string)
    title("figure \($(collect('a':'z')[fig_to_make])\)")
    if (fig_to_make != 3)
        xscale("log")
        #yscale("log")
    end
    yscale("log")
    legend(loc=2)
    savefig("julia_weissman_$(fig_to_make)_v2.pdf")
#    return low_sigma_results, high_sigma_results

#    return (sigmalist, UL_list, all_results)
end

#(low_sigma_results, high_sigma_results) = load_and_plot(1)
#load_and_plot(1)
load_and_plot(2)
load_and_plot(3)
load_and_plot(4)

# sigmalist, UL_list, nt_results = load_and_plot("neut_tunnel")
# UL_vals = unique(UL_list)
# #display(
# fig = figure()
#
# for (ii, i) in enumerate(UL_vals)
#     ax = errorbar(
#         sigmalist,
#         [mean(x) for x in nt_results[((ii-1)*length(sigmalist)+1):ii*length(sigmalist)]],
#         yerr = [1.0*std(x) for x in nt_results[((ii-1)*length(sigmalist)+1):ii*length(sigmalist)]],
#         label = "UL = $i"
#         )
# end
# ax = axes()
# legend(loc=2)
# ylabel(L"$\tau$")
# xlabel(L"$\sigma$")
# xscale("log")
# yscale("log")
# title("neutral tunneling")
# ylim([1e4,1e7])
# #savefig("neut_tunnel.pdf")
