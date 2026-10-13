"""Field certificate for the explicit degree-32 equation algebra.

This module uses only the displayed equations. It contains no curve search,
genus computation, divisor construction, or Riemann--Roch computation.
Run the full certificate with python3 verify.py from the repository root.
"""

F2 = GF(2)
K.<x> = FunctionField(F2)
R.<Y> = PolynomialRing(K)


def base_curve():
    equation = Y^4 + (x + 1) * Y^2 + (x^3 + x) * Y + x^7 + x^3
    return K.extension(equation, names="y")


# The explicit Artin-Schreier right-hand sides -- the same relations that
# check_F2_product_properties.sage uses -- as coefficient lists [y^0,y^1,y^2,y^3].
EXPLICIT_AS_RHS_COEFFS = [
    [
        "(x^9 + x^6 + x^5 + x^4 + x^3 + x^2)/(x^10 + x^4 + 1)",
        "(x^5 + x^4 + x^3 + x^2 + 1)/(x^10 + x^4 + 1)",
        "(x^4 + x^3 + x^2)/(x^10 + x^4 + 1)",
        "(x^3 + 1)/(x^10 + x^4 + 1)",
    ],
    [
        "(x^8 + x^5 + x^4 + x^3 + x^2 + x)/(x^10 + x^4 + 1)",
        "(x^6 + x^5 + x^2 + x)/(x^10 + x^4 + 1)",
        "(x^3 + x)/(x^10 + x^4 + 1)",
        "(x^4 + x^2 + x)/(x^11 + x^10 + x^5 + x^4 + x + 1)",
    ],
    [
        "(x^6 + x^4 + x^3 + x^2)/(x^10 + x^4 + 1)",
        "(x^5 + 1)/(x^10 + x^4 + 1)",
        "(x^3 + x^2 + x + 1)/(x^10 + x^4 + 1)",
        "(x^3 + x^2)/(x^10 + x^4 + 1)",
    ],
]


def certify_coordinate_algebra_is_a_field():
    """Certify that A_F = GF(2^128)(x)[y, z1, z2, z3]/(relations) is a field.

    This proves the field hypothesis for the equation-level norm bound.
    It also proves that GF(2) is the full constant field: the constant-field
    degree divides 32, and every nontrivial divisor of 32 has nontrivial gcd
    with 128, so a larger constant field would make A_F split. Over the perfect
    base field GF(2), full constants give geometric irreducibility. This is
    the point-count hypothesis; a singular rational point alone is insufficient.

    Method (nothing trusted beyond Sage's univariate irreducibility test over a
    rational function field, plus the explicit relations):
      1. build the explicit degree-32 algebra as a tower over GF(2);
      2. form the multiplication matrix of a primitive element t and take its
         characteristic polynomial m (degree 32 over GF(2)(x)); cross-check that
         m(t) = 0, so m genuinely presents the explicit relations;
      3. verify m is irreducible over the deployment field GF(2^128)(x).
    Step 3 is exactly the statement that A_F is a field.  (Sage cannot test
    irreducibility or genus directly on the relative tower, so we reduce to one
    irreducibility test over the rational function field, which it can do.)
    """
    L = base_curve()
    y = L.gen()
    M = L
    z_generators = []
    for coeffs in EXPLICIT_AS_RHS_COEFFS:
        x_in_M, y_in_M = M(x), M(y)
        rhs = sum(
            sage_eval(s, locals={"x": x_in_M}) * y_in_M^power
            for power, s in enumerate(coeffs)
        )
        Rz = PolynomialRing(M, names="Z")
        Z = Rz.gen()
        M = M.extension(Z^2 + Z - rhs)
        z_generators.append(M.gen())
    y_top = M(y)
    z1, z2, z3 = [M(g) for g in z_generators]

    def to_K_coords(element):
        # coordinates of a tower element in the 32 monomials y^j z1^a z2^b z3^c
        if element.parent() is K:
            return [K(element)]
        coords = []
        for c in element.list():
            coords += to_K_coords(c)
        return coords

    monomials = [
        y_top^j * z1^a * z2^b * z3^c
        for c in (0, 1) for b in (0, 1) for a in (0, 1) for j in range(4)
    ]
    basis_matrix = matrix(K, [to_K_coords(mono) for mono in monomials]).transpose()
    if not basis_matrix.is_invertible():
        raise RuntimeError("the 32 monomials are not a basis")

    t = y_top + z1 + M(x) * z2 + M(x)^2 * z3
    multiply_t = matrix(K, [to_K_coords(t * mono) for mono in monomials]).transpose()
    m = (basis_matrix.inverse() * multiply_t).charpoly()
    if m.degree() != 32:
        raise RuntimeError("multiplication charpoly has degree %s" % m.degree())

    # cross-check: m annihilates t, so m presents exactly the explicit algebra
    Rt = PolynomialRing(M, names="TT")
    TT = Rt.gen()
    if sum(M(c) * TT^i for i, c in enumerate(m.list()))(t) != 0:
        raise RuntimeError("characteristic polynomial does not annihilate t")
    if not m.is_irreducible():
        raise RuntimeError("t is not primitive / A is reducible over GF(2)(x)")

    # the certificate: A_F is a field iff m stays irreducible over GF(2^128)(x)
    F = GF(2^128)
    K_deploy = FunctionField(F, names="xd")
    xd = K_deploy.gen()

    def lift_to_deploy(coeff):
        num = sum(F(c) * xd^i for i, c in enumerate(coeff.numerator().list()))
        den = sum(F(c) * xd^i for i, c in enumerate(coeff.denominator().list()))
        return num / den

    Rd = PolynomialRing(K_deploy, names="Td")
    Td = Rd.gen()
    m_deploy = sum(lift_to_deploy(c) * Td^i for i, c in enumerate(m.list()))
    if not m_deploy.is_irreducible():
        raise RuntimeError(
            "A_F is NOT a field: minimal polynomial reducible over GF(2^128)(x)"
        )
    return ZZ(m.degree())
