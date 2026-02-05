# CPA self-consistency routine

"""
    cpa_loop(
        dispersion::Dispersion,
        ω::Float64 = 0.0;
        δ::Float64 = 0.1,
        x::Float64 = 0.0,
        ϵ::Float64 = 0.0,
        Σ::ComplexF64 = 0.0im,
        tol::Float64 = 0.1,
        maxiter::Int = 1000
    )

Calculate the CPA self-consistency loop for one frequency.
"""
function cpa_loop(
        dispersion::Dispersion,
        ω::Float64 = 0.0;
        δ::Float64 = 0.1,
        x::Float64 = 0.0,
        ϵ::Float64 = 0.0,
        Σ::ComplexF64 = 0.0im,
        tol::Float64 = 0.1,
        maxiter::Int = 1000
    )
    0.0 <= x <= 1.0 || throw(ArgumentError("x must be in [0, 1]"))
    δ > 0 || throw(ArgumentError("negative broadening"))
    tol > 0 || throw(ArgumentError("tol must be > 0"))
    maxiter >= 1 || throw(ArgumentError("maxiter must be >= 1"))
    G_loc = zero(ComplexF64)
    Σ = zero(ComplexF64)
    for _ in 1:maxiter
        # local GF
        G_loc = greens_function_local(dispersion, ω; Σ, δ)
        # take out self-energy
        𝒢_inv = inv(G_loc) + Σ
        𝒢 = inv(𝒢_inv)
        # resonant level model
        G_loc = (1 - x) * 𝒢 + x * inv(𝒢_inv - ϵ)
        # new self-energy
        Σ_new = 𝒢_inv - inv(G_loc)
        abs(Σ_new - Σ) < tol && break
        # new self-energy
        Σ = Σ_new
    end
    return G_loc, Σ
end

"""
    cpa_loop(
        dispersion::Dispersion,
        W::AbstractVector{Float64};
        δ::Float64 = 0.1,
        x::Float64 = 0.0,
        ϵ::Float64 = 0.0,
        Σ::Vector{ComplexF64} = zeros(ComplexF64, length(W)),
        tol::Float64 = 0.1,
        maxiter::Int = 1000
    )

Calculate the CPA self-consistency loop for given frequency grid.
"""
function cpa_loop(
        dispersion::Dispersion,
        W::AbstractVector{Float64};
        δ::Float64 = 0.1,
        x::Float64 = 0.0,
        ϵ::Float64 = 0.0,
        Σ::Vector{ComplexF64} = zeros(ComplexF64, length(W)),
        tol::Float64 = 0.1,
        maxiter::Int = 1000
    )
    G_loc = similar(W, ComplexF64)
    Σ_new = similar(W, ComplexF64)
    Threads.@threads for i in eachindex(W)
        G_loc[i], Σ_new[i] = cpa_loop(dispersion, W[i]; δ, x, ϵ, Σ = Σ[i], tol, maxiter)
    end
    return G_loc, Σ_new
end
