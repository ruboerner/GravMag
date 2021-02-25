"""
    dipole(m, r, rs)

Compute the magnetic anomaly of a magnetic dipole having a magnetic dipole moment `m` located at point `rs` observed in point `r`.
"""
function dipole(m::Array{T}, r::Array{T}, rs::Array{T}) where T<:Real
    R = norm(r - rs)
    return 4e-7 * pi / (4 * pi * R^3) *
        (
            3 * dot(r - rs, m) * (r - rs) / R^2 - m 
        )
end