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

w5 = CSV.read("weissman5a.csv")

fig = figure()
ax = []
for sigma in sort(unique(w5[:sigma]), rev=true)
    ws = w5[w5[:sigma].==sigma,:]
    ks = sort(unique(ws[:K]))
    push!(ax,
        errorbar(
                 [k for k in ks],
                 [mean(ws[ws[:K].==k,:crossing_time]) for k in ks],
                 yerr = vcat([quantile(ws[ws[:K].==k,:crossing_time], 0.25) for k in ks]',
                             [quantile(ws[ws[:K].==k,:crossing_time], 0.75) for k in ks]'),
                 alpha = 0.75,
                 label = L"$\sigma$ = " * string(sigma),
                 )
        )
end

analytical_prediction = [weissman5a_theory(2*k, [unique(w5[:mu1])[1],
                                               unique(w5[:mu2])[1],
                                               unique(w5[:s])[1],
                                               unique(w5[:delta])[1]]) for k in sort(unique(w5[:K]))];
plot(sort(unique(w5[:K])), analytical_prediction, c="k", ls="--")

xscale("log")
yscale("log")
ylabel(L"\tau")
display(fig)
#xlabel(indep_var_string)
