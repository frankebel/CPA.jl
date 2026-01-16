using CPA
using CairoMakie

function main()
    n = 20 # number of points

    ibz = ibz_2d(n)
    colorrange = extrema(ibz.multiplicity)

    f = Figure()
    ax = Axis(
        f[1, 1];
        title = L"multiplicity of $k$ points in the IBZ",
        xlabel = L"k_x",
        ylabel = L"k_y",
        # xticks = (0:(π / 4):π, [L"0", L"\frac{π}{4}", L"\frac{π}{2}", L"\frac{3π}{4}", L"π"]),
        xticks = (0:(π / 4):π, [L"0", L"π/4", L"π/2", L"3π/4", L"π"]),
        yticks = (0:(π / 4):π, [L"0", L"π/4", L"π/2", L"3π/4", L"π"]),
    )
    scatter!(ax, ibz.kgrid; color = ibz.multiplicity, colorrange = colorrange)
    Colorbar(f[1, 2], colorrange = colorrange, label = "multiplicity")

    save("ibz.pdf", f)

    return nothing
end

main()
