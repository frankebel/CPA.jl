# Green's function without impurities

using CairoMakie
using CPA


function main()
    W = range(-5; stop = 5, length = 500)
    δ = 0.01
    n = 1000
    ibz = ibz_2d(n)
    dispersion = dispersion_2d(ibz)
    g0 = greens_function_local(dispersion, W; δ)

    f = Figure()
    ax = Axis(
        f[1, 1];
        xlabel = L"ω/t",
        ylabel = L"A_0(\omega)",
    )
    lines!(ax, W, -imag(g0) ./ π)

    save("g0.pdf", f)

    return nothing
end

main()
