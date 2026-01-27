# Green's function

"""
    greens_function_local(
        dispersion::Dispersion,
        ω::Real;
        Σ::Number = 0.0im,
        δ::Real = 0.1,
    )

Calculate the local Green's function for one frequency.
"""
function greens_function_local(
        dispersion::Dispersion,
        ω::Real;
        Σ::Number = 0.0im,
        δ::Real = 0.1,
    )
    δ > 0 || throw(ArgumentError("negative broadening"))
    result = zero(ComplexF64)
    @fastmath @inbounds @simd for i in eachindex(dispersion.energy)
        G = inv(ω + im * δ - dispersion.energy[i] - Σ)
        result += dispersion.multiplicity[i] * G
    end
    result /= sum(dispersion.multiplicity)
    return result
end


"""
    greens_function_local(
        dispersion::Dispersion,
        W::Vector;
        Σ::AbstractVector = zeros(ComplexF64, length(W));
        δ::Real = 0.1,
    )

Calculate the local Green's function for a given frequency grid.

If no self-energy is given, it is assumed to be zero.
"""
function greens_function_local(
        dispersion::Dispersion,
        W::AbstractVector,
        Σ::AbstractVector = zeros(ComplexF64, length(W));
        δ::Real = 0.1,
    )
    length(W) == length(Σ) || throw(ArgumentError("length mismatch W, Σ"))
    result = Vector{ComplexF64}(undef, length(W))

    Threads.@threads for i in eachindex(W)
        @inbounds result[i] = greens_function_local(dispersion, W[i]; Σ = Σ[i], δ)
    end
    return result
end

function greens_function_local_resolvent(
        dispersion::Dispersion;
        Σ_H::Real = 0.0,
        Σ::PolesSum{Float64, Float64} = PolesSum([0.0], [0.0]),
    )
    Σ = remove_zero_weight(Σ)

    n = length(Σ) + 1
    n_tot = length(dispersion) * n
    loc_new = Vector{Float64}(undef, n_tot)
    wgt_new = Vector{Float64}(undef, n_tot)
    G_loc = PolesSum(loc_new, wgt_new)

    # Create tridiagonal matrix `T` using Householder transformations.
    # These do not touch the (1,1) element of the original matrix,
    # making it perfect for our use case.
    h = hessenberg!(Array(Σ)) # don't give symmetric information on purpose
    T = SymTridiagonal(diag(h.H), diag(h.H, -1)) # diagonal and first lower diagonal

    # Update the (1,1) element for each pole in Δ0 and diagonalize `T`.
    # Threads.@threads for i in eachindex(dispersion)
    for i in eachindex(dispersion)
        idx_low = 1 + n * (i - 1)
        idx_high = idx_low + n - 1
        bar = copy(T)
        bar[1, 1] = Σ_H + dispersion.energy[i]
        loc_new[idx_low:idx_high], U = eigen!(bar)
        U = h.Q * U # transform back
        wgt_new[idx_low:idx_high] = dispersion.multiplicity[i] .* abs2.(view(U, 1, :)) # multiply new weights with multiplicity
    end

    N_k = sum(dispersion.multiplicity)
    wgt_new .*= inv(N_k)
    sort!(G_loc)
    merge_degenerate_poles!(G_loc, eps())

    return G_loc
end
