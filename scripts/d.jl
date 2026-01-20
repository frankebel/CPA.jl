# exercise (d)

using CairoMakie
using CPA

function main()
    # parameters
    ω = 0.0 # only look at spectrum A(ω=0)
    δs = logrange(1.0e-2, 1.0e-8; length = 7)
    n = 1_500
    ibz = ibz_2d(n)
    dispersion = dispersion_2d(ibz)
    ϵs = range(3.8; stop = 4.2, length = 30)
    x = 0.5
    tol = 1.0e-10

    # plot A(ω=0) for varying ϵ, δ
    f = Figure()
    ax = Axis(
        f[1, 1];
        xlabel = L"\epsilon",
        ylabel = L"A(\omega=0)",
    )
    cmap = Colorbar(f[1, 2]; limits = extrema(δs), label = L"\delta", scale = log10)

    A0 = similar(ϵs)
    for (i, δ) in pairs(δs)
        @info "δ = $δ"
        for (j, ϵ) in pairs(ϵs)
            dispersion = dispersion_2d(ibz; μ = ϵ / 2) # shift chemical potential in each step
            G_loc, Σ = cpa_loop(dispersion, ω; δ, x, ϵ, tol)
            A0[j] = -imag(G_loc) / π
        end
        lines!(
            ax,
            ϵs,
            A0,
            color = log10(δ),
            colormap = :viridis,
            colorrange = log10.(extrema(δs)),
        )
    end
    save("d1.pdf", f)

    # plot A(ω=0, ϵ) for varying δ
    f = Figure()
    ax = Axis(
        f[1, 1];
        xlabel = L"\delta",
        ylabel = L"A(\omega=0)",
        xscale = log10,
        yscale = log10,
    )

    A0 = similar(δs)
    for ϵ in range(4.0; stop = 4.04, step = 0.01)
        @info "ϵ = $ϵ"
        for (i, δ) in pairs(δs)
            dispersion = dispersion_2d(ibz; μ = ϵ / 2) # shift chemical potential in each step
            G_loc, Σ = cpa_loop(dispersion, ω; δ, x, ϵ, tol)
            A0[i] = -imag(G_loc) / π
        end
        scatterlines!(
            ax,
            δs,
            A0,
            label = L"\epsilon=%$ϵ",
        )
    end
    axislegend(ax; position = :rb)
    save("d2.pdf", f)

    return nothing
end

main();
