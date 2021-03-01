function F_(u::T, v::T, w::T) where T<:Real
	R1 = sqrt(v^2 + w^2)
	R2 = sqrt(u^2 + w^2)
	R3 = sqrt(u^2 + v^2)
	R  = sqrt(u^2 + v^2 + w^2)
	return sign(u*v) * 
		(
			abs(u) * log( (abs(v) + R3) * R2 / (abs(u) * (abs(v) + R))) +
			abs(v) * log( (abs(u) + R3) * R1 / (abs(v) * (abs(u) + R))) +
			w * atan(abs(u * v) / w / R)
 		)
end

"""
   gbox(x, y, z, xextent, yextent, zextent)

Compute the vertical component of the gravitational attraction
rectangular box bounded by the axis-parallel planes defined 
by `xextent`, `yextent`, and `zextent`.
Observation point is given by the cartesian coordinates `x, y, z`.

Output is given in units of m/s^2.
"""
function gbox(x::T, y::T, z::T, x1::Array{T}, y1::Array{T}, z1::Array{T}) where T<:Real
	Vz = zero(eltype(T))
	for i in 1:2
        for j in 1:2
            for k in 1:2
                f = (-1.0)^(i+j+k)
                Vz += f * F_(x - x1[i], y - y1[j], z - z1[k])
            end
        end
    end
	return Vz * 6.6742e-11
end