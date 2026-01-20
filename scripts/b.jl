# exercise (b)

using CairoMakie
using CPA

function main()
    # parameters
    W = range(-7; stop = 7, length = 500)
    δ = 0.01
    n = 1000
    ibz = ibz_2d(n)
    ϵs = range(0; stop = 6, length = 20) # energies of disorder
    x = 0.5
    tol = 1.0e-6


    f = Figure()
    ax = Axis(
        f[1, 1];
        xlabel = L"\omega/t",
        ylabel = L"A(\omega)",
    )
    cmap = Colorbar(f[1, 2]; limits = extrema(ϵs), label = L"\epsilon")

    for (i, ϵ) in pairs(ϵs)
        dispersion = dispersion_2d(ibz; μ = ϵ / 2) # shift chemical potential in each step
        G_loc, Σ = cpa_loop(dispersion, W; δ, x, ϵ, tol)
        lines!(
            ax,
            W,
            -imag(G_loc) ./ π;
            color = ϵs[i],
            colormap = :viridis,
            colorrange = extrema(ϵs),
        )
    end

    save("b.pdf", f)
    return nothing
end

main();
