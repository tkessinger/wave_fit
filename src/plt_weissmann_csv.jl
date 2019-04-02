using CSV, PyPlot, Statistics

function weissman5a_theory(N, fixed_params)
    mu1, mu2, s, delta = fixed_params[1], fixed_params[2], fixed_params[3], fixed_params[4]
    #println("N = $N, mu1 = $mu1, mu2 = $mu2, s = $s, delta = $delta")
    if N < 1/sqrt(mu2*s)
        return 1.0/mu1 + 1.0/(N*mu2*(s+delta))
    elseif N < 1/mu1
        return 1.0/(N*mu1*sqrt(mu2*s))
    else
        return sqrt(pi/(2.0*N*mu1*mu2*s))
    end
end

function weissman5b_theory(N, fixed_params)
    mu1, mu2, s, delta = fixed_params[1], fixed_params[2], fixed_params[3], fixed_params[4]
    #println("N = $N, mu1 = $mu1, mu2 = $mu2, s = $s, delta = $delta")
    if N < 1/delta*log(1+(delta*(exp(delta)-1)/(mu2*s)))
        return 1.0/(N*mu1*(exp(delta)-1)/(exp(N*delta)-1)) + 1.0/(N*mu1*(s+delta))
        #return 1.0/(N*mu1*(exp(delta)-1)/(exp(N*delta)-1))
    elseif N < 2*delta^2/(pi*mu1*mu2*s)
        return delta/(N*mu1*mu2*s)
    else
        return sqrt(pi/(2.0*N*mu1*mu2*s))
    end
end

function plot_w5(data, theory, file = nothing)
    fig = figure()
    ax = []
    for sigma in sort(unique(data[:sigma]), rev=true)
        ds = data[data[:sigma].==sigma,:]
        ks = sort(unique(ds[:K]))
        push!(ax,
            errorbar(
                     [k for k in ks],
                     [mean(ds[ds[:K].==k,:crossing_time]) for k in ks],
                     yerr = vcat([quantile(ds[ds[:K].==k,:crossing_time], 0.25) for k in ks]',
                                 [quantile(ds[ds[:K].==k,:crossing_time], 0.75) for k in ks]'),
                     alpha = 0.75,
                     label = L"$\sigma$ = " * string(sigma),
                     )
            )
    end

    ks = sort(unique(data[:K]))
    mu1s = sort(unique(data[:mu1]))
    mu2s = sort(unique(data[:mu2]))
    ss  = sort(unique(data[:s]))
    deltas = sort(unique(data[:delta]))

    analytical_prediction = [theory(2*k, [mu1s[1], mu2s[1], ss[1], deltas[1]]) for k in ks];
    plot(sort(unique(data[:K])), analytical_prediction, c="k", ls="--")

    xscale("log")
    yscale("log")
    #xlabel(indep_var_string)
    ylabel(L"\tau")
    if file == nothing
        display(fig)
    end
end

w5a = CSV.read("weissman_fig5a_2019_03/weissman5a.csv")
w5b = CSV.read("weissman_fig5b_2019_03/weissman5b.csv")

plot_w5(w5a, weissman5a_theory)
plot_w5(w5b, weissman5b_theory)
