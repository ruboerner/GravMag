"""
    dipole(m, r, rs)

Compute the magnetic anomaly of a magnetic dipole having a magnetic dipole moment `m` located at point `rs` observed in point `r`.
"""
function dipole(m, r, rs)
    mu = pi * 4e-7
    R = norm(r - rs)
    B = mu / (4 * pi * R^3) *
        (
            3 * dot(r - rs, m) * (r - rs) / R^2 - m 
        )
    return B
end