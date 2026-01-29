module CPA

using LinearAlgebra
using RAS_DMFT

export
    # types
    # functions
    cpa_loop,
    cpa_loop_resolvent,
    dispersion_2d,
    greens_function_local,
    greens_function_local_resolvent,
    ibz_2d

include("brillouin_zone.jl")
include("dispersion.jl")
include("greens_function.jl")
include("loop.jl")

end
