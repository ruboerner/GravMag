module GravMag

using LinearAlgebra

# Type declarations
struct Point{T}
    x::T
    y::T
    z::T
end

abstract type Datum end

struct ScalarDatum <: Datum
    p::Point
    d::Float64
    datatype::String
    unit::String
end

ScalarDatum(p, d) = ScalarDatum(p, d, "dT", "nT")
struct VectorDatum <: Datum
	p::Point
	d::Array{Float64, 1}
	datatype::String
	unit::String
end

VectorDatum(p, d) = VectorDatum(p, d, "B", "nT")

include("mprism.jl")
include("dipole.jl")
include("tools.jl")

export Point, Datum, ScalarDatum, VectorDatum
export mprism, dipole, dircos


end # module
