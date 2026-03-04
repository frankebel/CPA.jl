# Quasiparticle weight
using CPA
using CSV
using CairoMakie
using DataFrames
using RAS_DMFT

ϵs = range(0.5; stop = 6, step = 0.5) # energies of disorder

Z = Float64[]
for ϵ in ϵs
    Σ = read_hdf5("self-energy_$(ϵ).h5", PolesSum{Float64, Float64})
    push!(Z, quasiparticle_weight(Σ, 1.0e-8))
end

f = Figure();
ax = Axis(
    f[1, 1];
    xlabel = L"\epsilon/t",
    ylabel = L"Z",
)
lines!(
    ax,
    ϵs,
    Z,
)
f
