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

struct VectorDatum <: Datum
	p::Point
	d::Array{Float64, 1}
	datatype::String
	unit::String
	function VectorDatum(p, d, datatype="B", unit="nT")
		new(p, d, datatype, unit)
	end
end

include("mprism.jl")
include("dipole.jl")
include("tools.jl")

export Datum, ScalarDatum, VectorDatum
export mprism, dipole, dircos


end # module
