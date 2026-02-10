# Green's function without impurities

# %%
# load modules
using CPA
using CairoMakie
using RAS_DMFT

# %%
# dispersion relation
grid = range(-5; stop = 5, length = 10_001)
n = 1000
ibz = ibz_2d(n)
dispersion = dispersion_2d(ibz)

# %%
# non-interacting GF
G0 = PolesSum(dispersion)
G0 = to_grid(G0, grid)
merge_small_weight!(G0, sqrt(eps()))

# %%
# plot it
W0 = range(-5; stop = 5, length = 2001)
W = range(-5; stop = 5, length = 2000) # without ω = 0
σ = 5.0e-3
δ = 5.0e-3
b = 0.4
f = Figure();
ax = Axis(
    f[1, 1];
    xlabel = L"ω/t",
    ylabel = L"A_0(\omega)",
    xticks = -4:2:4,
)
lines!(ax, W0, -imag(evaluate_gaussian(G0, W0, σ)) ./ π; label = L"Gaussian, $\sigma=%$σ$")
lines!(ax, W0, -imag(evaluate_lorentzian(G0, W0, σ)) ./ π; label = L"Lorentzian, $\delta=%$δ$")
lines!(ax, W, spectral_function_loggaussian(G0, W, b); label = L"log Gaussian, $b=%$b$")
axislegend(ax; position = :lt)

save("g0_resolvent.pdf", f)

return nothing
