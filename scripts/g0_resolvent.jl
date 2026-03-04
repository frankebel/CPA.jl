# Green's function without impurities

# %%
# load modules
using CPA
using CSV
using CairoMakie
using DataFrames
using RAS_DMFT

# %%
# dispersion relation
n = 10
ibz = ibz_2d(n)
dispersion = dispersion_2d(ibz)

# %%
# non-interacting GF
G0 = PolesSum(dispersion)
merge_small_weight!(G0, sqrt(eps()))
merge_degenerate_poles!(G0, 1.0e-4)

# %%
# plot it
W = range(-5; stop = 5, length = 2001)
σ = δ = 2.0e-1
f = Figure();
ax = Axis(
    f[1, 1];
    xlabel = L"ω/t",
    ylabel = L"residues of $G_0$",
    xticks = -4:2:4,
)
plot!(ax, locations(G0), weights(G0))
ax = Axis(
    f[1, 2];
    xlabel = L"ω/t",
    ylabel = L"A_0(\omega)",
    xticks = -4:2:4,
)
lines!(ax, W0, -imag(evaluate_gaussian(G0, W0, σ)) ./ π; label = L"Gaussian, $\sigma=%$σ$")
lines!(ax, W0, -imag(evaluate_lorentzian(G0, W0, σ)) ./ π; label = L"Lorentzian, $\delta=%$δ$")
axislegend(ax; position = :lt)
f

# %%
# export
df = DataFrame()
df[!, "locations"] = locations(G0)
df[!, "weights"] = weights(G0)
CSV.write("g0.csv", df)

df = DataFrame()
df[!, "W"] = W
df[!, "gaussian"] = -imag(evaluate_gaussian(G0, W, σ)) ./ π
df[!, "lorentzian"] = -imag(evaluate_lorentzian(G0, W, δ)) ./ π
CSV.write("g0_broadened.csv", df)
