# methods related to irreducible Brillouin zone IBZ

"""
    struct IBZ{D}

Represent the irreducible Brillouin zone (IBZ).
"""
struct IBZ{D}
    kgrid::Vector{NTuple{D, Float64}}
    multiplicity::Vector{Int}

    function IBZ{D}(kgrid2d, multiplicity) where {D}
        length(kgrid2d) == length(multiplicity) ||
            throw(ArgumentError("mismatch list of momenta and their associated multiplicity"))
        return new{D}(kgrid2d, multiplicity)
    end
end

"""
    IBZ(kgrid::Vector{<:NTuple{D, Float64}}, multiplicity) where {D}

return a new `IBZ` instance  with given momentum grid `kgrid`
with each point having multiplicity `multiplicity`)
"""
function IBZ(kgrid::Vector{<:NTuple{D, Float64}}, multiplicity) where {D}
    return IBZ{D}(kgrid, multiplicity)
end

"""
    ibz_2d(n::Int, a::Real = 1)

Return a irreducible Brillouin zone on a 2d square lattice with ``n`` steps
in the interval ``[0, π/a]``.
"""
function ibz_2d(n::Int, a::Real = 1)
    kgrid = kgrid_IBZ_2d(n, a)
    multiplicity = multiplicity_IBZ_2d(kgrid, last(last(kgrid)))
    return IBZ(kgrid, multiplicity)
end


"""
    klist(n::Int, a::Real = 1)

Create an equidistant grid from `0` to `π/a` containing `n` points.
"""
function klist(n::Int, a::Real = 1)
    n >= 1 || throw(ArgumentError("must have at least 1 point"))
    return range(0; length = n, stop = π / a)
end

"""
    kgrid_IBZ_2d(n::Int, a::Real = 1)

Return a vector containing all pairs ``(k_x, k_y)`` in the IBZ.
For each dimension create an equidistant grid from `0` to `π/a` containing `n` points.
"""
function kgrid_IBZ_2d(n::Int, a::Real = 1)
    ks = klist(n, a)
    result = Vector{NTuple{2, eltype(ks)}}(undef, n * (n + 1) ÷ 2) # Gauss sum for number of points
    for x in eachindex(ks)
        for y in 1:x
            result[(x - 1) * x ÷ 2 + y] = (ks[x], ks[y])
        end
    end
    return result
end

"""
    multiplicity_IBZ_2d(kgrid2d::AbstractVector, klast::Real)

Find the multiplicity of each k-point in a 2d square lattice.
The variable `klast` is the highest value a momentum component can have.
"""
function multiplicity_IBZ_2d(kgrid::AbstractVector, klast::Real)
    #    M
    #   /|
    #  / |
    # Γ--X
    result = Vector{Int}(undef, length(kgrid))
    for (i, kpoint) in pairs(kgrid)
        if iszero(kpoint[1]) && iszero(kpoint[2])
            # Γ
            multiplicity = 1
        elseif kpoint[1] == klast && iszero(kpoint[2])
            # X
            multiplicity = 2
        elseif kpoint[1] == kpoint[2] == klast
            # M
            multiplicity = 1
        elseif iszero(kpoint[2])
            # Γ -- X
            multiplicity = 4
        elseif kpoint[1] == klast
            # X -- M
            multiplicity = 4
        elseif kpoint[1] == kpoint[2]
            # Γ -- M
            multiplicity = 4
        else
            # nothing special
            multiplicity = 8
        end
        result[i] = multiplicity
    end
    return result
end
