module CPA

using LinearAlgebra
using RAS_DMFT

export
    # types
    # functions
    cpa_loop,
    dispersion_2d,
    greens_function_local,
    ibz_2d

include("brillouin_zone.jl")
include("dispersion.jl")
include("greens_function.jl")
include("loop.jl")

end
