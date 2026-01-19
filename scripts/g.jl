# Green's function

using CairoMakie
using CPA

function main()
    # parameters
    W = range(-5; stop = 10, length = 500)
    δ = 0.01
    n = 1000
    ibz = ibz_2d(n)
    dispersion = dispersion_2d(ibz)

    x = 0.5
    ϵ = 4.0
    tol = 1.0e-6

    # G_loc, Σ = cpa_loop(dispersion, 0.0; δ, x, ϵ, tol) # single frequency
    G_loc, Σ = cpa_loop(dispersion, W; δ, x, ϵ, tol) # frequency grid

    f = Figure()
    ax = Axis(
        f[1, 1];
        xlabel = L"ω/t",
        ylabel = L"A(\omega)",
    )
    lines!(ax, W, -imag(G_loc) ./ π)

    save("g.pdf", f)

    return nothing
end

main()
