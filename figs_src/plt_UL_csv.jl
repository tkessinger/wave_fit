using CSV, PyPlot, Statistics

ct = CSV.read("output/parallel_sims/var_UL_weissman5b.csv")
ct = [ct;CSV.read("output/parallel_sims/var_UL_weissman5b_2.csv")]

#sigma_tags = reverse([L"$5 \times 10^{-2}$",L"$1 \times 10^{-2}$",L"$5 \times 10^{-8}$"])
#sigma_tags = reverse([L"$5 \times 10^{-1}$", L"$1 \times 10^{-1}$", L"$5 \times 10^{-2}$", L"$1 \times 10^{-2}$",L"$5 \times 10^{-8}$", L"$1 \times 10^{-8}$"])

sigma_tags = Dict(
5e-1 => L"$5 \times 10^{-1}$",
1e-1 => L"$1 \times 10^{-1}$",
5e-2 => L"$5 \times 10^{-2}$",
1e-2 => L"$1 \times 10^{-2}$",
5e-8 => L"$5 \times 10^{-8}$",
1e-8 => L"$1 \times 10^{-8}$"
)

for (si, sigma) in enumerate(sort(unique(ct[:sigma])))
    fig = figure()
    println(sigma)
    ax = []
    for (ULi, UL) in enumerate(sort(unique(ct[:UL])))
        cts = ct[(ct[:sigma].==sigma) .& (ct[:UL].==UL),:]
        ks = sort(unique(cts[:K]))
        push!(ax,
            errorbar(
                     ks,
                     [mean(cts[cts[:K].==k,:crossing_time]) for k in ks],
                     yerr = vcat([quantile(cts[cts[:K].==k,:crossing_time], 0.25) for k in ks]',
                                 [quantile(cts[cts[:K].==k,:crossing_time], 0.75) for k in ks]'),
                     alpha = 0.75,
                     label = L"$UL = $" * string(UL)
                     )
            )
    end
    title(L"$\sigma = $" * sigma_tags[sigma])
    xscale("log")
    yscale("log")
    ylabel(L"\tau")
    xlabel(L"N")
    legend(loc=2)
    tight_layout()
    display(fig)
    savefig("var_UL_$sigma.pdf")
end

#xlabel(indep_var_string)
