# Resolvent with linear grid

# %%
# load modules
using CPA
using CairoMakie
using RAS_DMFT

# %%
# dispersion relation
grid = range(-10; stop = 10, length = 1001)
n = 1000
ϵ = 1.5
μ = ϵ / 2
x = 0.5

# non-interacting GF
ibz = ibz_2d(n)
dispersion = dispersion_2d(ibz; μ)
G0 = PolesSum(dispersion)
G0 = to_grid(G0, grid)
merge_small_weight!(G0, sqrt(eps()))

# %%
# CPA
maxiter = 15
G, Σ_H, Σ = cpa_loop_resolvent(G0, grid; x, ϵ, maxiter)
write_hdf5("G_$(ϵ).h5", G0)
write_hdf5("self-energy_Hartree_$(ϵ).h5", Σ_H)
write_hdf5("self-energy_$(ϵ).h5", Σ)

# %%
# plot it
W = range(-10; stop = 10, length = 2001);
σ = δ = 5.0e-2;
f = Figure(; size = (1600 / 2, 900 / 2));
ax = Axis(
    f[1, 1];
    xlabel = L"ω/t",
    ylabel = L"A_\mathrm{loc}(\omega)",
);
lines!(ax, W, -imag(evaluate_gaussian(G, W, σ)) ./ π);

ax = Axis(
    f[1, 2];
    xlabel = L"ω/t",
    ylabel = L"-\mathrm{Im}~\Sigma(\omega)",
    #=xticks = -4:2:4,=#
);
lines!(ax, W, -imag(evaluate_gaussian(Σ, W, σ)));
display(f)

# %%
# moments
moment(G)
moment(Σ)
ϵ^2 / 4
