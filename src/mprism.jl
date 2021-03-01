"""
    mprism(x, y, z, xextent, yextent, zextent, M)

Compute the magnetic anomaly of a prism bounded by the axis-parallel planes defined
by `xextent`, `yextent`, and `zextent`.
Observation point is given by cartesian coordinates `x, y, z`.
The magnetization of the body is defined by the three components of `M` measured in units of A/m.
"""
function mprism(x::T, y::T, z::T, x1::Array{T}, y1::Array{T}, z1::Array{T}, M::Array{T}) where T<:Real
    dZ = zero(eltype(T))
    dX = zero(eltype(T))
    dY = zero(eltype(T))
    f = zero(eltype(T))
    for i in 1:2
        for j in 1:2
            for k in 1:2
                f = (-1.0)^(i+j+k)
                dX += f * (
                M[1] * g2(y - y1[j], z - z1[k], x - x1[i]) +
                M[2] * g1(x - x1[i], z - z1[k], y - y1[j]) +
                M[3] * g1(x - x1[i], y - y1[j], z - z1[k])
                )

                dY += f * (
                M[1] * g1(x - x1[i], z - z1[k], y - y1[j]) +
                M[2] * g2(x - x1[i], z - z1[k], y - y1[j]) +
                M[3] * g1(y - y1[j], x - x1[i], z - z1[k])
                )

                dZ += f * (
                M[1] * g1(x - x1[i], y - y1[j], z - z1[k]) +
                M[2] * g1(y - y1[j], x - x1[i], z - z1[k]) +
                M[3] * g2(x - x1[i], y - y1[j], z - z1[k]))
            end
        end
    end
    dX *= 1e-7
    dY *= 1e-7
    dZ *= 1e-7
    return [dX, dY, dZ]
end

function g1(u::T, v::T, w::T) where T<:Real
    return sign(v) * log(
        sqrt(u^2 + w^2) /
        (abs(v) + sqrt(u^2 + v^2 + w^2))
        )
end

function g2(u::T, v::T, w::T) where T<:Real
    return atan(u * v / (w * sqrt(u^2 + v^2 + w^2)))
end

"""
    mprism_(x, y, z, xextent, yextent, zextent, M)

Compute the magnetic anomaly B of a prism bounded by the axis-parallel planes defined
by `xextent`, `yextent`, and `zextent`.
Observation point is given by cartesian coordinates `x, y, z`.
The magnetization of the body is defined by the three components of `M` measured in units of A/m.
"""
function mprism_(x::T, y::T, z::T, x1::Array{T}, y1::Array{T}, z1::Array{T}, M::Array{T}) where T<:Real
	G11 = G2(y, z, x, y1, z1, x1)
	G12 = G1(x, z, y, x1, z1, y1)
	G13 = G1(x, y, z, x1, y1, z1)
	#G21 = G12
	G22 = G2(x, z, y, x1, z1, y1)
	G23 = G1(y, x, z, y1, x1, z1)
	#G31 = G13
	#G32 = G23
	G33 = G2(x, y, z, x1, y1, z1)
	
	G = [G11 G12 G13
		G12 G22 G23
		G13 G23 G33]
		
	B = G * M 
end

"""
    mprism!(B, x, y, z, xextent, yextent, zextent, M)

Compute the magnetic anomaly B of a prism bounded by the axis-parallel planes defined
by `xextent`, `yextent`, and `zextent`.
Observation point is given by cartesian coordinates `x, y, z`.
The magnetization of the body is defined by the three components of `M` measured in units of A/m.

Note that this function overrides the value of the argument `B`.
"""
function mprism!(B::Array{T}, x::T, y::T, z::T, x1::Array{T}, y1::Array{T}, z1::Array{T}, M::Array{T}) where T<:Real
	G11 = G2(y, z, x, y1, z1, x1)
	G12 = G1(x, z, y, x1, z1, y1)
	G13 = G1(x, y, z, x1, y1, z1)
	#G21 = G12
	G22 = G2(x, z, y, x1, z1, y1)
	G23 = G1(y, x, z, y1, x1, z1)
	#G31 = G13
	#G32 = G23
	G33 = G2(x, y, z, x1, y1, z1)
	
	B .= [G11 G12 G13
		G12 G22 G23
		G13 G23 G33] * M
	return nothing	
	# B .= G * M 
end

function G2(x::T, y::T, z::T, x1::Array{T, 1}, y1::Array{T, 1}, z1::Array{T, 1}) where T<:Real
        return -g2(y - y1[1], z - z1[1], x - x1[1]) +
         g2(y - y1[1], z - z1[2], x - x1[1]) +
         g2(y - y1[2], z - z1[1], x - x1[1]) -
         g2(y - y1[2], z - z1[2], x - x1[1]) +
         g2(y - y1[1], z - z1[1], x - x1[2]) -
         g2(y - y1[1], z - z1[2], x - x1[2]) -
         g2(y - y1[2], z - z1[1], x - x1[2]) +
         g2(y - y1[2], z - z1[2], x - x1[2])
end


function G1(x::T, y::T, z::T, x1::Array{T, 1}, y1::Array{T, 1}, z1::Array{T, 1}) where T<:Real
        return -g1(y - y1[1], z - z1[1], x - x1[1]) +
         g1(y - y1[1], z - z1[2], x - x1[1]) +
         g1(y - y1[2], z - z1[1], x - x1[1]) -
         g1(y - y1[2], z - z1[2], x - x1[1]) +
         g1(y - y1[1], z - z1[1], x - x1[2]) -
         g1(y - y1[1], z - z1[2], x - x1[2]) -
         g1(y - y1[2], z - z1[1], x - x1[2]) +
         g1(y - y1[2], z - z1[2], x - x1[2])
end
