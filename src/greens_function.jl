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

    # Create tridiagonal matrix `T` using Householder transformations.
    # These do not touch the (1,1) element of the original matrix,
    # making it perfect for our use case.
    h = hessenberg!(Array(Σ)) # don't give symmetric information on purpose
    T = SymTridiagonal(diag(h.H), diag(h.H, -1)) # diagonal and first lower diagonal

    # chunks for multithreading
    chunk_size = max(1, length(dispersion) ÷ Threads.nthreads())
    chunk_data = Iterators.partition(eachindex(dispersion.energy), chunk_size)

    tasks = map(chunk_data) do chunk
        Threads.@spawn begin
            # allocate once
            bar = copy(T)
            baz = copy(h.Q)

            # return data
            loc = Float64[]
            wgt = Float64[]

            for i in chunk
                copyto!(bar, T)
                # Update the (1,1) element for each pole and diagonalize `T`.
                bar[1, 1] = Σ_H + dispersion.energy[i]
                Λ, U = eigen!(bar)
                mul!(baz, h.Q, U) # transform back

                # write into result
                append!(loc, Λ)
                v = baz[1, :]
                @. v = abs2(v) * dispersion.multiplicity[i] # new weights scaled by original
                append!(wgt, v)
            end

            return loc, wgt
        end
    end

    states = fetch.(tasks)

    # merge to long vectors
    loc_new = mapreduce(i -> i[1], vcat, states)
    wgt_new = mapreduce(i -> i[2], vcat, states)
    G_loc = PolesSum(loc_new, wgt_new)

    # normalization
    N_k = sum(dispersion.multiplicity)
    wgt_new .*= inv(N_k)
    sort!(G_loc)
    merge_degenerate_poles!(G_loc, eps())

    return G_loc
end

"""
    PolesSum(dispersion::Dispersion)

Create a `PolesSum` instance from a given dispersion relation `dispersion`.
"""
function RAS_DMFT.PolesSum(dispersion::Dispersion)
    loc = Vector{Float64}(undef, length(dispersion))
    wgt = Vector{Float64}(undef, length(dispersion))
    G = PolesSum(loc, wgt)

    @inbounds for i in eachindex(dispersion)
        loc[i] = dispersion.energy[i]
        wgt[i] = dispersion.multiplicity[i]
    end

    # normalize
    N_k = sum(dispersion.multiplicity)
    wgt .*= inv(N_k)

    # cleanup
    sort!(G)
    merge_degenerate_poles!(G, 10 * eps())

    return G
end
