using CSV, PyPlot, Statistics
#
# function weissman5a_theory(N, fixed_params)
#     mu1, mu2, s, delta = fixed_params[1], fixed_params[2], fixed_params[3], fixed_params[4]
#     #println("N = $N, mu1 = $mu1, mu2 = $mu2, s = $s, delta = $delta")
#     if N < 1/sqrt(mu2*s)
#         return 1.0/mu1 + 1.0/(N*mu2*(s+delta))
#     elseif N < 1/mu1
#         return 1.0/(N*mu1*sqrt(mu2*s))
#     else
#         return sqrt(pi/(2.0*N*mu1*mu2*s))
#     end
# end
#
# w5 = CSV.read("weissman5a.csv")
#
# fig = figure()
# ax = []
# for sigma in sort(unique(w5[:sigma]), rev=true)
#     ws = w5[w5[:sigma].==sigma,:]
#     ks = sort(unique(ws[:K]))
#     push!(ax,
#         errorbar(
#                  [k for k in ks],
#                  [mean(ws[ws[:K].==k,:crossing_time]) for k in ks],
#                  yerr = vcat([quantile(ws[ws[:K].==k,:crossing_time], 0.25) for k in ks]',
#                              [quantile(ws[ws[:K].==k,:crossing_time], 0.75) for k in ks]'),
#                  alpha = 0.75,
#                  label = L"$\sigma$ = " * string(sigma),
#                  )
#         )
# end
#
# analytical_prediction = [weissman5a_theory(2*k, [unique(w5[:mu1])[1],
#                                                unique(w5[:mu2])[1],
#                                                unique(w5[:s])[1],
#                                                unique(w5[:delta])[1]]) for k in sort(unique(w5[:K]))];
# plot(sort(unique(w5[:K])), analytical_prediction, c="k", ls="--")
#
# xscale("log")
# yscale("log")
# ylabel(L"\tau")
# display(fig)
#xlabel(indep_var_string)

using CSV, PyPlot, Statistics, Base.MathConstants, PyCall

rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.size"] = 14
rcParams["lines.linewidth"] = 2

function seqfix(N, fixed_params)
    mu1, mu2, s, delta = fixed_params[1], fixed_params[2], fixed_params[3], fixed_params[4]

    return 1.0/(N*mu1*(exp(delta)-1)/(exp(N*delta)-1)) + 1.0/(N*mu2*(s+delta))
end

function w5_theory(N, fixed_params)
    mu1, mu2, s, delta = fixed_params[1], fixed_params[2], fixed_params[3], fixed_params[4]

    # eq 25 and 22 from Weissman et al (2009)
    p = (-delta + sqrt(delta^2 + 4*mu2*s))/2
    tau = log(2/(1+delta/sqrt(delta^2+4*mu2*s)))/p

    if delta < 0 && - N * delta > 1
        # beneficial single mutant (also from tunnel eq 27)
        return 1.0/(N*mu1*p) + tau + eulergamma/s
    elseif delta < 2 * sqrt(mu2*s)
        # neutral single mutant
        if N < 1/sqrt(mu2*s)
            # Sequential fixation
            return 1.0/(N*mu1*(exp(delta)-1)/(exp(N*delta)-1)) + 1.0/(N*mu2*(s+delta))
        elseif N < 1/mu1
            # Tunneling (eq 27)
            return 1.0/(N*mu1*p) + tau + eulergamma/s
        elseif N < 2*s/(pi*mu1*mu2)
            # Semi-determinstic
            return sqrt(pi/(2.0*N*mu1*mu2*s))
        else
            # Deterministic
            return log((s+delta)/(N*mu1*mu2))
        end
    else
        # deleterious single mutant
        if N < 1/delta*log(1+delta*(exp(delta)-1)/(mu2*s))
            # Sequential fixation
            return 1.0/(N*mu1*(exp(delta)-1)/(exp(N*delta)-1)) + 1.0/(N*mu2*(s+delta))
        elseif delta < s
            if N < 2*delta^2/(pi*mu1*mu2*s)
                # Tunneling (eq 27)
                return 1.0/(N*mu1*p) + tau + eulergamma/s
            elseif N < 2*s/(pi*mu1*mu2)
                # Semi-determinstic
                return sqrt(pi/(2.0*N*mu1*mu2*s))
            else
                # Deterministic
                return log((s+delta)/(N*mu1*mu2))
            end
        else
            if N < delta/(mu1*mu2)
                # Tunneling (eq 27)
                return 1.0/(N*mu1*p) + tau + eulergamma/s
            else
                # Deterministic
                return log((s+delta)/(N*mu1*mu2))
            end
        end
    end
end

function plot_w5(data, theory, var; figlabel = "", xlog = true, ylog = true, ymax = nothing,
                 loc = 3, xtic = nothing, ytic = nothing, file = nothing, show_legend = false)
    fig = figure()
    ax = []
    #sigma_tags = [L"0.5", L"0.1", L"0.05", L"0.01", L"10^{-8}"]
    sigma_tags = [L"5 \times 10^{-1}", L"1 \times 10^{-1}", L"5 \times 10^{-2}", L"1 \times 10^{-2}", L"1 \times 10^{-8}"]
    for (si, sigma) in enumerate(sort(unique(data[:sigma]), rev=true))
        ds = data[data[:sigma].==sigma,:]
        vs = sort(unique(ds[var]))
        push!(ax,
            errorbar(
                     [v for v in vs],
                     [mean(ds[ds[var].==v,:crossing_time]) for v in vs],
                     yerr = vcat([quantile(ds[ds[var].==v,:crossing_time], 0.25) for v in vs]',
                                 [quantile(ds[ds[var].==v,:crossing_time], 0.75) for v in vs]'),
                     alpha = 0.75,
                     label = L"$\sigma$ = " * sigma_tags[si]
                     )
            )
    end

    ks = sort(unique(data[:K]))
    mu1s = sort(unique(data[:mu1]))
    mu2s = sort(unique(data[:mu2]))
    ss  = sort(unique(data[:s]))
    deltas = sort(unique(data[:delta]))
    ULs = sort(unique(data[:UL]))
    sigmas = sort(unique(data[:sigma]))

    if var == :K
        indep_var_string = L"N"
        analytical_prediction_drift = [theory(2*k, [mu1s[1], mu2s[1], ss[1], deltas[1]]) for k in ks];
        plot(sort(unique(data[var])), analytical_prediction_drift, c="k", ls="--")
        # t2 = (N,σ) -> sqrt(24*log(N*σ))/(N*σ)*N
        # analytical_prediction_draft = [seqfix(2*t2(k,1), [mu1s[1], mu2s[1], ss[1], 1e-4]) for k in ks];
        # plot(sort(unique(data[var])), analytical_prediction_draft, c="k", ls="-.")
    elseif var == :delta
        indep_var_string = L"\delta"
        analytical_prediction_drift = [theory(2*ks[1], [mu1s[1], mu2s[1], ss[1], d]) for d in deltas];
        plot(sort(unique(data[var])), analytical_prediction_drift, c="k", ls="--")
    elseif var == :mu2
        indep_var_string = L"\mu_2"
        analytical_prediction_drift = [theory(2*ks[1], [mu1s[1], mu2, ss[1], deltas[1]]) for mu2 in mu2s];
        plot(sort(unique(data[var])), analytical_prediction_drift, c="k", ls="--")
    end

    # draftf = k -> (6*log(k*(ULs[1]*sigmas[3]^2/2)^(1/3))/(ULs[1]*sigmas[3]^2))^(1/3)
    # analytical_prediction_draft = [theory(2*draftf(k), [mu1s[1], mu2s[1], ss[1], deltas[1]]) for k in ks];
    # plot(sort(unique(data[:K])), analytical_prediction_draft, c="k", ls="-.")

    if xlog
        xscale("log")
    end
    if ylog
        yscale("log")
    end

    if xtic != nothing
        xticks(xtic)
    end

    if ytic != nothing
        yticks(ytic)
    end

    if ylim != nothing
        ylim(0, ymax)
    end

    fig.text(0.07,.95, figlabel, size=18)#, weight="bold")

    xlabel(indep_var_string)
    ylabel(L"\tau")
    if show_legend == true
        legend(loc=loc, fontsize=12)
    end
    tight_layout()
    if file == nothing
        display(fig)
    else
        savefig(file)
    end
end

w5a = CSV.read("5a_2019_03/weissman5a.csv")
w5b = CSV.read("5b_2019_03/weissman5b.csv")
w5c = CSV.read("5c_2019_04/weissman5c.csv")
w5d = CSV.read("5d_2019_04/weissman5d.csv")

plot_w5(w5a, w5_theory, :K, file="julia_weissman_5a.pdf", figlabel="(a)", show_legend = true)
plot_w5(w5b, w5_theory, :K, figlabel = "(b)", file="julia_weissman_5b.pdf")
plot_w5(w5c, w5_theory, :delta, figlabel = "(c)", xlog=false, loc=4, xtic=[-0.01, -0.005, 0.0, 0.005, 0.01], file="julia_weissman_5c.pdf")
plot_w5(w5d, w5_theory, :mu2, figlabel = "(d)", file="julia_weissman_5d.pdf")
