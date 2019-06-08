using CSV, PyPlot, Statistics, IterTools

#load the crossing times
#ct = CSV.read("output/parallel_sims/grid.csv")
#ct = CSV.read("output/parallel_sims/grid_sims_s01_UL1_highmu.csv")
ct = CSV.read("output/parallel_sims/grid_sims_s01_UL1.csv")

#load the sweep times for comparison
st = CSV.read("output/parallel_sims/grid_sims_sweeps_UL1_ver2.csv")

#define the tick labels and sigma titles
sigma_tags = [L"$5 \times 10^{-8}$",L"$5 \times 10^{-2}$"]
#delta_tags = [L"$10^{"* string(x) * L"}$" for x in range(-6.0,stop=-2.0,step=0.5)]
#N_tags = reverse([L"$10^{"* string(x) * L"}$" for x in range(3,stop=6,step=1.0/3)])
delta_tags = [L"$10^{"* string(round(Int64,x)) * L"}$" for x in range(-6.0,stop=-2.0,step=0.5)]
N_tags = reverse([L"$10^{"* string(round(Int64,x)) * L"}$" for x in range(3,stop=6,step=1.0/3)])

ks = sort(unique(ct[:K])) #get unique values of K (i.e., N)
deltas = sort(unique(ct[:delta])) #get unique values of δ
valley_s = sort(unique(ct[:s]))[1]
sigma_vals = sort(unique(ct[:sigma])) #get unique values of σ
ks_and_deltas = collect(Iterators.product(ks,deltas))
sweep_s = -1.0*sort(unique(st[:delta]))
ks_and_svals = collect(Iterators.product(ks,sweep_s))

ratios = []

tauplots = false
alphaplots = false
ratioplots = true

for (si, sigma) in enumerate(sigma_vals)
    cts = ct[ct[:sigma].==sigma,:] #just pick the ct values for the current sigma value

    #average all the crossing time values for fixed k and δ
    mean_cts = [mean(cts[(cts[:K].==k) .& (cts[:delta].==delta),:crossing_time]) for (k, delta) in ks_and_deltas]
    for (k, delta) in (ks_and_deltas)
        println("$k, $delta, $(mean(cts[(cts[:K].==k) .& (cts[:delta].==delta),:crossing_time])), " *
        "$(length(cts[(cts[:K].==k) .& (cts[:delta].==delta),:crossing_time]))")
    end
    mean_cts = reverse(mean_cts, dims=1) #so that N increases along the y axis
    if tauplots
        fig, ax1 = subplots(1,1)
        imshow(max.(0,log10.(mean_cts)),vmin=0,vmax=ceil(maximum(log10.(mean_cts)))+1)
        cbar = colorbar()

        #labels are chosen so as to avoid ugly fractional exponents
        ax1.set_xticks(0:2:8)
        ax1.set_xticklabels(delta_tags[1:2:9])
        ax1.set_xlabel(L"$\delta$")
        ax1.set_yticks(0:3:9)
        ax1.set_yticklabels(N_tags[1:3:10])
        ax1.set_ylabel(L"N")
        cbar.set_label(L"$\log_{10} \tau$")
        #title(L"$\sigma = 10^{"*log10(sigma)"}$")
        title(L"$\sigma = $" * sigma_tags[si])
        tight_layout()

        show(fig)
    end

    sts = st[st[:sigma].==sigma,:]
    for (svi, sval) in enumerate(sweep_s)
        #mean_sts = [mean(sts[(sts[:K].==k) .& (sts[:delta].==-1.0*s),:crossing_time]) for (k, s) in ks_and_svals]
        mean_sts = [mean(sts[(sts[:K].==k) .& (sts[:delta].==-1.0*sval), :crossing_time]) for (k, delta) in ks_and_deltas]
        mean_sts = reverse(mean_sts, dims=1)
        if alphaplots
            fig, ax1 = subplots(1,1)
            imshow(log10.(mean_cts./mean_sts),vmin=-2,vmax=2,cmap="RdBu")
            #imshow(max.(0,log10.(mean_sts)),vmin=0,vmax=ceil(maximum(log10.(mean_cts)))+1)

            cbar = colorbar()

            ax1.set_xticks(0:2:8)
            ax1.set_xticklabels(delta_tags[1:2:9])
            ax1.set_xlabel(L"$\delta$")
            ax1.set_yticks(0:3:9)
            ax1.set_yticklabels(N_tags[1:3:10])
            ax1.set_ylabel(L"N")
            cbar.set_label(L"$\log_{10} \alpha = \tau_{\mathrm{valley}}/\tau_{\mathrm{sweep}}$")
            #title(L"$\sigma = 10^{"*log10(sigma)"}$")
            title(L"$\sigma = $" * sigma_tags[si] * L"$, s_{\mathrm{valley}} = $" * string(valley_s)
                 * L"$, s_{\mathrm{sweep}} = $" * string(sval))
            #tight_layout()
            show(fig)
        end
        if sval == 0.001
            push!(ratios, mean_cts./mean_sts)
        end
    end

end

if ratioplots
    fig, ax1 = subplots(1,1)
    imshow(log10.(ratios[2]./ratios[1]), vmin=-1, vmax = 1, cmap="PuOr")
    cbar = colorbar()
    ax1.set_xticks(0:2:8)
    ax1.set_xticklabels(delta_tags[1:2:9])
    ax1.set_xlabel(L"$\delta$")
    ax1.set_yticks(0:3:9)
    ax1.set_yticklabels(N_tags[1:3:10])
    ax1.set_ylabel(L"N")
    cbar.set_label(L"$\log_{10} (\alpha_{\sigma_2}/\alpha_{\sigma_1})$")
    #title(L"$\sigma = 10^{"*log10(sigma)"}$")
    title("ratio of "* L"$\alpha$" *" values, " * L"$s_{\mathrm{sweep}} = 0.01$")
    #tight_layout()
    show(fig)
end

#xscale("log")
#yscale("log")
#ylabel(L"\tau")
#show(fig)
#xlabel(indep_var_string)
