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
    length(W) == length(Σ) || throw(DimensionMismatch("length mismatch between W, Σ"))
    G_loc = similar(W, ComplexF64)
    Σ_new = similar(W, ComplexF64)
    Threads.@threads for i in eachindex(W)
        G_loc[i], Σ_new[i] = cpa_loop(dispersion, W[i]; δ, x, ϵ, Σ = Σ[i], tol, maxiter)
    end
    return G_loc, Σ_new
end

"""
    cpa_loop_resolvent(
        dispersion::Dispersion,
        grid::AbstractVector{Float64};
        x::Real = 0.0,
        ϵ::Real = 0.0,
        Σ_H::Real = 0.0,
        Σ::PolesSum{Float64, Float64} = PolesSum([0.0], [0.0]),
        maxiter::Int = 100,
    )

Calculate the CPA loop in the resolvent formalism.
"""
function cpa_loop_resolvent(
        G0::PolesSum,
        grid::AbstractVector{Float64};
        x::Real = 0.0,
        ϵ::Real = 0.0,
        Σ_H::Real = 0.0,
        Σ::PolesSum{Float64, Float64} = PolesSum([0.0], [0.0]),
        maxiter::Int = 100,
    )
    # local GF
    G_loc = greens_function_local_resolvent(G0; Σ_H, Σ)
    G_loc = to_grid(G_loc, grid)
    @info "length G_loc $(length(G_loc))"
    for it in 1:maxiter
        @info "iteration $it"
        # local GF
        G_loc = greens_function_local_resolvent(G0; Σ_H, Σ)
        G_loc = to_grid(G_loc, grid)
        # impurity GF
        a_0, G_loc_inv = inv(G_loc)
        a_0 -= Σ_H
        𝒢0_inv = G_loc_inv - Σ
        merge_negative_weight!(𝒢0_inv)
        merge_small_weight!(𝒢0_inv, eps())
        # 𝒢0
        foo = Array(𝒢0_inv)
        foo[1, 1] = a_0
        F = eigen!(foo)
        wgt = F.vectors[1, :]
        map!(abs2, wgt)
        𝒢0 = PolesSum(F.values, wgt)
        # resonant level model (RLM)
        # G_loc = (1 - x) * 𝒢0 + x * (𝒢0^{-1} - ϵ)^{-1}
        #       = P1           + P2
        P1 = copy(𝒢0)
        weights(P1) .*= 1 - x
        foo = Array(𝒢0_inv)
        foo[1, 1] = a_0 + ϵ
        F = eigen!(foo)
        wgt = F.vectors[1, :]
        map!(abs2, wgt)
        wgt .*= x
        P2 = PolesSum(F.values, wgt)
        G_loc = P1 + P2
        # new self-energy
        a_loc, G_loc_inv = inv(G_loc)
        Σ_H = a_loc - a_0
        Σ = G_loc_inv - 𝒢0_inv
        merge_negative_weight!(Σ)
        merge_small_weight!(Σ, sqrt(eps()))
    end

    return G_loc, Σ_H, Σ
end
