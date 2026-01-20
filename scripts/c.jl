# exercise (c)

using CairoMakie
using CPA

function main()
    # parameters
    W = range(-5; stop = 5, length = 500)
    δ = 0.01
    n = 1000
    ibz = ibz_2d(n)
    dispersion = dispersion_2d(ibz)
    ϵ = 0.001
    x = 0.01
    tol = 1.0e-6

    G_loc, Σ = cpa_loop(dispersion, W; δ, x, ϵ, tol)

    f = Figure(; size = (900, 300))
    ax = Axis(f[1, 1]; xlabel = L"\omega/t", ylabel = L"\mathrm{Re}~\Sigma(\omega) - x\epsilon")
    refline = x * ϵ
    lines!(ax, W, real(Σ) .- refline)
    ax = Axis(f[1, 2]; xlabel = L"\omega/t", ylabel = L"-\mathrm{Im}~\Sigma(\omega)")
    lines!(ax, W, -imag(Σ); label = L"-\mathrm{Im}~\Sigma(\omega)")
    lines!(ax, W, -imag(greens_function_local(dispersion, W; δ)) .* x .* ϵ^2; linestyle = :dash, label = L"\pi x\epsilon^2 A_0(\omega)")
    axislegend(ax; position = :lt)

    save("c.pdf", f)
    return nothing
end

main();
