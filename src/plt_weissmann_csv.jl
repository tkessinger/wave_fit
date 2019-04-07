using CSV, PyPlot, Statistics

function w5_theory(N, fixed_params)
    mu1, mu2, s, delta = fixed_params[1], fixed_params[2], fixed_params[3], fixed_params[4]

    # eq 25 from Weissman et al (2009)
    p = (-delta + sqrt(delta^2 + 4*mu2*s))/2

    if delta < 0 && - N * delta > 1
        # beneficial single mutant (also from tunnel eq 27)
        return 1.0/(N*mu1*p)
    elseif delta < 2 * sqrt(mu2*s)
        # neutral single mutant
        if N < 1/sqrt(mu2*s)
            # Sequential fixation
            return 1.0/(N*mu1*(exp(delta)-1)/(exp(N*delta)-1)) + 1.0/(N*mu2*(s+delta))
        elseif N < 1/mu1
            # Tunneling (eq 27)
            return 1.0/(N*mu1*p)
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
                return 1.0/(N*mu1*p)
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
                return 1.0/(N*mu1*p)
            else
                # Deterministic
                return log((s+delta)/(N*mu1*mu2))
            end
        end
    end
end

function plot_w5(data, theory, var; xlog = true, ylog = true, ymax = nothing, loc = 3, file = nothing)
    fig = figure()
    ax = []
    for sigma in sort(unique(data[:sigma]), rev=true)
        ds = data[data[:sigma].==sigma,:]
        vs = sort(unique(ds[var]))
        push!(ax,
            errorbar(
                     [v for v in vs],
                     [mean(ds[ds[var].==v,:crossing_time]) for v in vs],
                     yerr = vcat([quantile(ds[ds[var].==v,:crossing_time], 0.25) for v in vs]',
                                 [quantile(ds[ds[var].==v,:crossing_time], 0.75) for v in vs]'),
                     alpha = 0.75,
                     label = L"$\sigma$ = " * string(sigma)
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
        analytical_prediction_drift = [theory(2*k, [mu1s[1], mu2s[1], ss[1], deltas[1]]) for k in ks];
        plot(sort(unique(data[var])), analytical_prediction_drift, c="k", ls="--")
    elseif var == :delta
        analytical_prediction_drift = [theory(2*ks[1], [mu1s[1], mu2s[1], ss[1], d]) for d in deltas];
        plot(sort(unique(data[var])), analytical_prediction_drift, c="k", ls="--")
    elseif var == :mu2
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

    if ylim != nothing
        ylim(0, ymax)
    end

    #xlabel(indep_var_string)
    ylabel(L"\tau")
    legend(loc=loc)
    if file == nothing
        display(fig)
    end
end

w5a = CSV.read("weissman_fig5a_2019_03/weissman5a.csv")
w5b = CSV.read("weissman_fig5b_2019_03/weissman5b.csv")
w5c = CSV.read("weissman_fig5c_2019_04/weissman5c.csv")
w5d = CSV.read("weissman_fig5d_2019_04/weissman5d.csv")

plot_w5(w5a, w5_theory, :K)
plot_w5(w5b, w5_theory, :K)
plot_w5(w5c, w5_theory, :delta, xlog=false, loc=4)
plot_w5(w5d, w5_theory, :mu2)
