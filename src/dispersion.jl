# dispersion relation

"""
    struct Dispersion{R <: Real}

Represent the dispersion relation in the irreducible Brillouin zone (IBZ) together
with its multiplicity.
"""
struct Dispersion
    energy::Vector{Float64}
    multiplicity::Vector{Int}

    function Dispersion(energy, multiplicity)
        length(energy) == length(multiplicity) ||
            throw(ArgumentError("mismatch list of energies and their associated multiplicity"))
        return new(energy, multiplicity)
    end
end

"""
    dispersion_2d(ibz::IBZ; μ::Real = 0, t::Real = 1)

Calculate the dispersion relation

```math
ϵ_k = -2t \\left( \\cos(k_x) + \\cos(k_y) \\right)
```

for the given irreducible Brillouin zone `ibz`.
"""
function dispersion_2d(ibz::IBZ; μ::Real = 0, t::Real = 1)
    energy = map(k -> energy_2d(k, t), ibz.kgrid)
    energy .-= μ
    return Dispersion(energy, ibz.multiplicity)
end

energy_2d(k::NTuple{2, <:Real}, t::Real = 1) = -2 * t * (cos(k[1]) + cos(k[2]))

Base.eachindex(dispersion::Dispersion) = eachindex(dispersion.energy)

Base.length(dispersion::Dispersion) = length(dispersion.energy)
