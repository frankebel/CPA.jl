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
