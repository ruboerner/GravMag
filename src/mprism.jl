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
