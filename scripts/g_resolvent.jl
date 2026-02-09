# One CPA step

# %%
# load modules
using CPA
using CairoMakie
using RAS_DMFT

# %%
# dispersion relation
grid = range(-10; stop = 10, length = 101)
n = 100
ϵ = 8.0
μ = ϵ / 2
x = 0.5

ibz = ibz_2d(n)
dispersion = dispersion_2d(ibz; μ)

# %%
# CPA

maxiter = 100
G_loc1, Σ_H1, Σ1 = cpa_loop_resolvent(dispersion, grid; x, ϵ, maxiter)
G_loc2, Σ2 = cpa_loop(dispersion, grid; x, ϵ, maxiter);

# plot it
W = range(-10; stop = 10, length = 2001);
σ = δ = 5.0e-2;
f = Figure(; size = (1600, 900));
ax = Axis(
    f[1, 1];
    xlabel = L"ω/t",
    ylabel = L"A_\mathrm{loc}(\omega)",
);
lines!(ax, W, -imag(evaluate_lorentzian(G_loc1, W, δ)) ./ π; label = "resolvent");
lines!(ax, grid, -imag(G_loc2) ./ π; label = "finite broadening");
axislegend(ax; position = :lt);

ax = Axis(
    f[1, 2];
    xlabel = L"ω/t",
    ylabel = L"\mathrm{Im}~\Sigma(\omega)",
    #=xticks = -4:2:4,=#
);
lines!(ax, W, -imag(evaluate_lorentzian(Σ1, W, δ)) ./ π; label = "resolvent");
lines!(ax, grid, -imag(Σ2) ./ π; label = "finite broadening");
axislegend(ax; position = :lt);

display(f)

# %%
# moments

moment(G_loc1)
moment(G_loc2, grid)
moment(Σ1)
moment(Σ2, grid)
