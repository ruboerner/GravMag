function dircos(incl, decl)
    xincl = incl * pi / 180.0
    xdecl = decl * pi / 180.0
    a = cos(xincl) * cos(xdecl)
    b = cos(xincl) * sin(xdecl)
    c = sin(xincl)
    return a, b, c
end