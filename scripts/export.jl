# export data
using CPA
using CSV
using CairoMakie
using DataFrames
using RAS_DMFT

ϵs = range(0; stop = 6, step = 0.5) # energies of disorder
W = range(-7; stop = 7, length = 501)
σ = 3.0e-2

# %%
# Green's function
df = DataFrame()
df[!, "W"] = W

f = Figure();
ax = Axis(
    f[1, 1];
    xlabel = L"\omega/t",
    ylabel = L"A(\omega)",
)
cmap = Colorbar(f[1, 2]; limits = extrema(ϵs), label = L"\epsilon")

for (i, ϵ) in pairs(ϵs)
    G = read_hdf5("G_$(ϵ).h5", PolesSum{Float64, Float64})
    g = -imag(evaluate_gaussian(G, W, σ)) ./ π
    df[!, "$ϵ"] = g
    lines!(
        ax,
        W,
        g,
        color = ϵs[i],
        colormap = :viridis,
        colorrange = extrema(ϵs),
    )
end
f
# CSV.write("g_linear.csv", df)

# %%
# self-energy
df = DataFrame()
df[!, "W"] = W
f = Figure();
ax = Axis(
    f[1, 1];
    xlabel = L"\omega/t",
    ylabel = L"-\mathrm{Im}\Sigma",
)
cmap = Colorbar(f[1, 2]; limits = extrema(ϵs), label = L"\epsilon")

for (i, ϵ) in pairs(ϵs)
    iszero(ϵ) && continue
    Σ = read_hdf5("self-energy_$(ϵ).h5", PolesSum{Float64, Float64})
    s = -imag(evaluate_gaussian(Σ, W, σ))
    df[!, "$ϵ"] = s
    lines!(
        ax,
        W,
        s,
        color = ϵs[i],
        colormap = :viridis,
        colorrange = extrema(ϵs),
    )
end
xlims!(-4, 4)
f
CSV.write("self-energy_linear.csv", df)
