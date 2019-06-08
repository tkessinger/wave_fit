using CSV, PyPlot, Statistics

ct = CSV.read("output/parallel_sims/var_sigma_delta.csv")
#ct = [ct;CSV.read("output/parallel_sims/var_sigma_delta_2.csv")]


#sigma_tags = [L"$5 \times 10^{-8}$",L"$5 \times 10^{-2}$"]

fig = figure()
ax = []
for sval in sort(unique(ct[:s]))
    for (si, sigma) in enumerate(sort(unique(ct[:sigma])))
        cts = ct[(ct[:sigma].==sigma) .& (ct[:s].==sval),:]
        deltas = sort(unique(cts[:delta]))
        push!(ax,
            errorbar(
                     [delta for delta in deltas],
                     [mean(cts[cts[:delta].==delta,:crossing_time]) for delta in deltas],
                     yerr = vcat([quantile(cts[cts[:delta].==delta,:crossing_time], 0.25) for delta in deltas]',
                                 [quantile(cts[cts[:delta].==delta,:crossing_time], 0.75) for delta in deltas]'),
                     alpha = 0.75,
                     label = L"$\sigma = $" * sigma_tags[si] * L"$, s = $" * string(sval)
                     )
            )
    end
end

xscale("log")
yscale("log")
ylabel(L"\tau")
xlabel(L"\delta")
legend(loc=2)
tight_layout()
display(fig)
savefig("julia_sigma_delta.pdf")
#xlabel(indep_var_string)
