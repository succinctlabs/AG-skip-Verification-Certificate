"""Deterministic certificate for the paper's equation-level claims.

Run through ../verify.py. No search outputs or precomputed success flags are used.
The mathematical implications behind the numerical checks are in README.md.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import time

import sage.all as sage
from sage.repl.preparse import preparse
from sage.version import version as sage_version

AUDIT = Path(__file__).resolve().parent
START = time.monotonic()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def report(name, value):
    print(f"{name}: {value}", flush=True)


def load_namespace(name):
    path = AUDIT / name
    namespace = dict(vars(sage), __name__="certificate_module", __file__=str(path))
    exec(compile(preparse(path.read_text()), str(path), "exec"), namespace)
    return namespace


def check_paper_data(a, paper=None):
    snapshot = json.loads((AUDIT / "paper_instance.json").read_text())
    rows = [sum(int(v) << j for j, v in enumerate(row)) for row in a["extension_matrix_R"].rows()]
    require(rows == [int(v, 16) for v in snapshot["R_rows_u64"]], "R differs from the printed paper matrix")
    eta = a["K"](a["product_denominator"])
    require(eta == sum(a["x"]**i for i in snapshot["eta_nonzero_powers"]), "eta differs from the paper")
    if paper:
        source = Path(paper)
        if source.is_dir():
            source /= "explicit_instance.tex"
        text = source.read_text()
        printed = re.findall(r"(\d+)\s*&\s*\\texttt\{(0x[0-9a-fA-F]+)\}", text)
        require([int(i) for i, _ in printed] == list(range(222)), "paper must contain 222 indexed matrix rows")
        require([int(v, 16) for _, v in printed] == rows, "supplied paper matrix differs from R")
        equation = text.split(r"\label{eq:eta}", 1)[1].split(r"\end{equation}", 1)[0]
        expression = equation.split("=", 1)[1].replace("\n", "").replace(" ", "").rstrip(".")
        # Accept only the displayed sum of x powers, so an unexpected term fails.
        terms = expression.split("+")
        powers = []
        for term in terms:
            match = re.fullmatch(r"x(?:\^\{(\d+)\}|\^(\d+))?", term)
            require(match is not None, f"unrecognized eta term: {term}")
            powers.append(int(match.group(1) or match.group(2) or 1))
        require(eta == sum(a["x"]**i for i in powers), "supplied paper eta differs from the evaluator")
        report("supplied_paper_matrix_and_eta", "ok")
    report("paper_matrix_rows", len(rows))
    report("paper_eta", "ok")
    return snapshot


def check_structure(a):
    F2 = sage.GF(2)
    R = a["extension_matrix_R"]
    require(R.base_ring() == F2 and R.dimensions() == (222, 64), "wrong extension matrix type")
    require(R[:64, :] == sage.identity_matrix(F2, 64), "R is not systematic")
    expected = [("finite", i, j) for j in range(3) for i in range(64)]
    expected += [("finite", i, 3) for i in range(30)]
    require(a["product_coordinate_profile"] == expected, "wrong 64+64+64+30 coordinate order")
    require(a["cover_basis"] == [(u, v, w) for u in (0, 1) for v in (0, 1) for w in (0, 1)], "wrong cover monomial order")
    # Verify B using independent formal inputs: coefficient comparison covers
    # all 222^2 ordered basis pairs without repeatedly evaluating zero terms.
    polynomial_ring = sage.PolynomialRing(F2, names=[f"u{i}" for i in range(222)] + [f"v{i}" for i in range(222)])
    variables = polynomial_ring.gens()
    left = sage.vector(polynomial_ring, variables[:222])
    right = sage.vector(polynomial_ring, variables[222:])
    expected_output = sage.vector(polynomial_ring, 222)
    offsets = [0, 64, 128, 192]
    for point in range(64):
        count = 4 if point < 30 else 3
        for left_order in range(count):
            for right_order in range(count - left_order):
                expected_output[offsets[left_order + right_order] + point] += left[offsets[left_order] + point] * right[offsets[right_order] + point]
    require(a["product_rule"](left, right) == expected_output, "symbolic block convolution mismatch")
    report("systematic_extension", "ok")
    report("product_systematicity_and_binary_preservation", "ok")
    report("block_convolution_symbolic_check", "ok (all 222^2 coefficients)")


def degree_at_infinity(f):
    return f.numerator().degree() - f.denominator().degree()


def check_bounds(a):
    x, y, K = a["x"], a["y"], a["K"]
    nu = 32
    d = (x + 1) * (x**10 + x**4 + 1)
    eta = K(a["product_denominator"])
    require(d.numerator().degree() == 11 and eta.numerator().degree() == 49, "unexpected denominator degrees")
    # Monic triangular relations give a free rank-32 module after inverting d.
    degrees = [7]
    for rhs in a["artin_schreier_rhs"]:
        coefficients = rhs.list()
        denominator = sage.lcm([f.denominator() for f in coefficients])
        require((d / denominator).denominator() == 1, "structural denominator is not inverted by d")
        degrees.append(max([denominator.degree() + 2] + [(f * denominator).numerator().degree() + j for j, f in enumerate(coefficients) if f]))
    require(sorted(degrees) == [7, 12, 12, 13], "cleared equation degrees differ from the paper")
    delta = sage.prod(degrees)

    # y has pole weight <=7/4: every lower term in its monic relation has
    # weight <=7. The AS right sides have nonpositive weight, so z_i have
    # weight <=0. This holds at every place above infinity.
    wy = sage.QQ(7) / 4
    base_coefficients = a["base_curve_polynomial"].list()
    require(all(not f or degree_at_infinity(f) + j * wy <= 4 * wy for j, f in enumerate(base_coefficients[:-1])), "invalid y pole weight")
    for rhs in a["artin_schreier_rhs"]:
        require(all(not f or degree_at_infinity(f) + j * wy <= 0 for j, f in enumerate(rhs.list())), "invalid z pole weight")

    h = x**5 + x**2 + 1
    require(h**2 == x**10 + x**4 + 1, "wrong square denominator")
    scalings = [h, (x + 1)*h, h]
    for scale, rhs in zip(scalings, a["artin_schreier_rhs"]):
        # (scale*z)^2 + scale*(scale*z) = scale^2*rhs is monic and integral.
        require(all(f.denominator() == 1 for f in (scale**2 * rhs).list()), "finite-pole clearing failed")
    clearing = sage.prod(scalings)
    require(clearing.numerator().degree() == 16, "finite-pole clearing degree is not 16")

    entries = [a["polynomial_function_coordinates"](gamma) for gamma in a["gamma_numerators"]]
    require(len(entries) == 222 and all(len(row) == 32 for row in entries), "wrong E array dimensions")
    require(all(f.denominator() == 1 for row in entries for f in row), "E contains nonpolynomial coefficients")
    max_x = max(f.numerator().degree() for row in entries for f in row if f)
    weighted = max(f.numerator().degree() + (j % 4) * wy for row in entries for j, f in enumerate(row) if f)
    require(max_x == 49, "E maximum x degree differs from the paper")
    # The paper's conservative bound uses degree<=49 and y-exponent<=3,
    # even though inspection of the exact support gives weighted degree 49.
    coarse_weight = max_x + 3 * wy
    norm_bound = nu * (clearing.numerator().degree() + coarse_weight)
    require(norm_bound == 2248, "norm bound differs from the paper")

    # Raw determinant bound: cleared entries have weighted degree <=80.
    # Individual entries can have degree >80; row/column weights cancel in
    # every determinant term. Check the reduction table itself, not metadata.
    raw_clearing = (x**10 + x**4 + 1)**3 * (x + 1)
    require((raw_clearing / clearing).denominator() == 1, "raw denominator does not clear finite poles")
    maxima = [max((f.numerator().degree() if f else -1) for f in column) for column in zip(*entries)]
    raw_entry_degree = 0
    weighted_entry_degree = 0
    for i, coefficient_degree in enumerate(maxima):
        if coefficient_degree < 0:
            continue
        for j in range(32):
            cover_product = a["reduce_monomial_product"](a["cover_basis"][i // 4], a["cover_basis"][j // 4])
            product = [c * y**(i % 4 + j % 4) for c in cover_product]
            for row, coefficient in enumerate(a["polynomial_function_coordinates"](product)):
                cleared = coefficient * raw_clearing
                require(cleared.denominator() == 1, "raw matrix denominator not cleared")
                if cleared:
                    degree = coefficient_degree + cleared.numerator().degree()
                    raw_entry_degree = max(raw_entry_degree, degree)
                    weighted_entry_degree = max(weighted_entry_degree, degree + (row % 4 - j % 4)*wy)
    require(weighted_entry_degree <= 80, "raw weighted determinant bound failed")
    raw_norm = nu * 80
    raw_zeros = nu * raw_norm + nu * (11 + 49)
    require(raw_zeros == 83840, "raw zero bound differs from the paper")

    q = sage.ZZ(2)**128
    require(q > 4*delta**2, "Lang-Weil size hypothesis failed")
    # Round Delta^(1/3) UP using integers; no floating-point bound decisions.
    cube_root_ceiling = delta.nth_root(3, truncate_mode=True)[0]
    if cube_root_ceiling**3 < delta:
        cube_root_ceiling += 1
    error_upper = (delta-1)*(delta-2)*2**64 + 5*delta**4*cube_root_ceiling + 11*delta + 49*nu
    require(error_upper < 2**92, "point-count error exceeds the paper's bound")
    n_lower = q - 2**92
    require(n_lower > norm_bound, "injectivity bound failed")
    require(raw_zeros * 2**111 < n_lower, "111-bit bound failed")
    require(norm_bound * 2**116 < n_lower, "116-bit bound failed")
    values = dict(field_degree=nu, equation_degrees=[int(v) for v in degrees], delta=int(delta),
                  denominator_degree=49, finite_pole_clearing_degree=16,
                  evaluator_max_x_degree=int(max_x), evaluator_weighted_degree=str(weighted),
                  raw_matrix_entry_degree_bound=int(raw_entry_degree),
                  raw_matrix_weighted_degree_bound=str(weighted_entry_degree),
                  raw_norm_bound=int(raw_norm), raw_zero_bound=int(raw_zeros),
                  norm_bound=int(norm_bound), point_count_lower_bound=int(n_lower),
                  equation_bound_bits=116, geometric_119_bit_claim="not checked by this certificate")
    for key, value in values.items():
        report(key, value)
    return values


def check_direct_base_evaluator(a):
    # Direct basis functions have a different common denominator. Compare
    # rational functions exactly, then check that the shortcut is defined on Y.
    base_denominator = a["K"](a["base_denominator"])
    eta = a["K"](a["product_denominator"])
    require(base_denominator != 0 and base_denominator.denominator() == 1, "invalid direct base denominator")
    require(eta.denominator() == 1, "invalid product denominator")
    for factor, _ in base_denominator.numerator().factor():
        require(eta.numerator() % factor == 0, "direct base denominator has a pole on Y")
    require(len(a["base_gamma_numerators"]) == 64, "wrong direct base evaluator dimension")
    for i, direct in enumerate(a["base_gamma_numerators"]):
        require(len(direct) == 8, "wrong direct base monomial count")
        via_product = a["combine_product_evaluator_numerators"](a["extension_matrix_R"].column(i))
        require(all(eta * direct[j] == base_denominator * via_product[j] for j in range(8)), f"direct base evaluator differs at basis vector {i}")
    report("direct_base_evaluator_symbolic_checks", 64)


def smoke_test(a):
    # Deterministic PRNG seed for repeatability. The certificate does not rely
    # on these evaluations. Exercise F2 and nonbinary GF(2^128) messages.
    sage.set_random_seed(20261006)
    field = a["RandomField"]
    for binary in (True, False):
        for _ in range(3):
            point = a["random_curve_point"]()
            x0, y0 = point[:2]
            require(y0**4 + (x0+1)*y0**2 + (x0**3+x0)*y0 + x0**7+x0**3 == 0, "sampler base equation failed")
            for z, rhs in zip(point[2:], a["artin_schreier_rhs"]):
                require(z**2+z == a["eval_base_function"](rhs, point), "sampler cover equation failed")
            base = sage.GF(2) if binary else field
            left = sage.vector(base, [base.random_element() for _ in range(64)])
            right = sage.vector(base, [base.random_element() for _ in range(64)])
            a["check_product_identity_at"](point, left, right)
            require(a["base_code_evaluate_direct"](left, point) == a["base_code_evaluate"](left, point), "direct base runtime check failed")
    report("binary_and_extension_field_smoke_checks", 6)


def input_hashes(paper=None):
    paths = [AUDIT.parent / "verify.py", AUDIT / "certify.py", AUDIT / "paper_instance.json"]
    paths += sorted(AUDIT.glob("*.sage"))
    result = {str(p.relative_to(AUDIT.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    if paper:
        source = Path(paper)
        if source.is_dir():
            source /= "explicit_instance.tex"
        result["supplied_paper"] = hashlib.sha256(source.read_bytes()).hexdigest()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper")
    parser.add_argument("--smoke", action="store_true", help="skip the field and 4096-product certificates")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    hashes = input_hashes(args.paper)
    report("SageMath", sage_version)
    report("mode", "SMOKE ONLY" if args.smoke else "FULL DETERMINISTIC CERTIFICATE")
    a = load_namespace("check_F2_product_properties.sage")
    snapshot = check_paper_data(a, args.paper)
    check_structure(a)
    rank_R, rank_E = a["check_nondegeneracy"]()
    report("R_rank", rank_R)
    report("product_evaluator_numerator_rank", rank_E)
    check_direct_base_evaluator(a)
    bounds = check_bounds(a)
    if not args.smoke:
        curve = load_namespace("check_F2_coordinate_field.sage")
        require(curve["EXPLICIT_AS_RHS_COEFFS"] == a["artin_schreier_rhs_coeffs"], "field and product checks use different curves")
        require(curve["base_curve"]().polynomial().list() == a["base_curve_polynomial"].list(), "base curve equations differ")
        report("coordinate_algebra_check", "starting")
        degree = curve["certify_coordinate_algebra_is_a_field"]()
        require(degree == 32, "wrong coordinate field degree")
        # Its constant-field degree c divides 32. For every possible c>1,
        # GF(2^c) and GF(2^128) have nontrivial intersection; the tensor
        # product would split. Deployment-field irreducibility excludes all c.
        require(all(sage.gcd(c, 128) > 1 for c in sage.divisors(degree) if c > 1), "constant-field inference does not apply")
        report("coordinate_algebra_degree", degree)
        report("full_constant_field_GF2", "ok")
        checked = a["check_symbolic_basis_products"]()
        require(checked == 4096, "incomplete product certificate")
        report("symbolic_basis_product_checks", checked)
    smoke_test(a)
    require(input_hashes(args.paper) == hashes, "certificate inputs changed during the run; rerun on a stable checkout")
    result = {"status": "smoke_only" if args.smoke else "certified", "sage_version": sage_version,
              "paper_revision": snapshot["paper_revision"], "bounds": bounds,
              "elapsed_seconds": round(time.monotonic()-START, 2),
              "inputs_sha256": hashes}
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    report("result", "SMOKE ONLY; field and product certificates skipped" if args.smoke else "PASS: equation-level certificate, 116-bit bound")
    report("elapsed_seconds", result["elapsed_seconds"])


if __name__ == "__main__":
    main()
