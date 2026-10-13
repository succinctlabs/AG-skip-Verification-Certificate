"""
Map-level product-property audit for the F2 genus-95 curve-code artifacts.

This script checks only the abstract maps that the code needs:

    base message --R--> product-code coordinates --evaluate at P--> field value

Base evaluation is defined as the product-code evaluator composed with R.
The optional direct base evaluator is checked against that definition; the
complete certificate also checks all 64 basis functions symbolically.

It is organized for human audit:

    1. Common utilities shared by both checks.
    2. Nondegeneracy/rank checks for R and the exported evaluator functions.
    3. A symbolic check proving the product identity exactly.
    4. Random checks exercising the same maps at GF(2^128)-points.
    5. Runtime cross-checks of the direct base evaluator.

The complete distance certificate is run by python3 verify.py; the proof and
claim-to-check mapping are in F2_human_audit/README.md.

Run from the project root with:

    $ sage F2_human_audit/check_F2_product_properties.sage

or from this folder with:

    $ sage check_F2_product_properties.sage

For a fast smoke test that skips the exact symbolic audit:

    $ sage F2_human_audit/check_F2_product_properties.sage --quick

Quick mode still runs the nondegeneracy/rank checks.

The new sage.cli launcher (conda-forge / pip Sage) does not pass script flags
through; under it, use environment variables instead:

    $ AUDIT_QUICK=1 sage F2_human_audit/check_F2_product_properties.sage

(also AUDIT_SYMBOLIC_ONLY=1 and AUDIT_RANDOM_TRIALS=N).

Sage may create a preparsed `.sage.py` file next to this script before running.
That file is generated cache material and is ignored by git.

Notation used below:

    gamma_i
        The numerator function for the i-th product-code coordinate.  The file
        F2_genus95_candidate_product_evaluator_polynomials.sage stores these as
        gamma_numerator_coeffs, reduced in the eight monomials z1^a z2^b z3^c.

    H
        The common denominator product_denominator exported by the evaluator.

    Gamma_w
        For a product-code coordinate vector w, the linear combination

            Gamma_w = sum_i w_i * gamma_i.

        Thus

            product_code_evaluate(w, P) = Gamma_w(P) / H(P)

        wherever the exported rational functions are defined.

For a base basis vector e_i, the symbolic check writes G_i for Gamma_{R*e_i}.
For the product-rule vector product_rule(R*e_i, R*e_j), it writes G_ij for the
corresponding Gamma.  The checked identity is

    G_i * G_j = H * G_ij

after reducing both sides by the three Artin-Schreier cover equations.
"""


# -----------------------------------------------------------------------------
# Common utilities
# -----------------------------------------------------------------------------

import argparse
import os
import sys
import time

EXTENSION_DEGREE_OVER_F2 = 128
RANDOM_TRIALS = 20
MAX_POINT_ATTEMPTS = 20000
SYMBOLIC_PROGRESS_INTERVAL = 512

def _locate_audit_dir():
    """Directory holding the sibling .sage files this script loads.

    The classic `sage` launcher sets __file__ to this script, so
    dirname(__file__) is correct.  The new sage.cli launcher (conda-forge / pip
    Sage) execs a temporary preparsed copy, so __file__ points into a temp
    directory; fall back to the working directory in that case.
    """
    sentinel = "R_matrix_F2_genus95_candidate.sage"
    candidates = []
    try:
        candidates.append(os.path.dirname(os.path.abspath(__file__)))
    except NameError:
        pass
    cwd = os.getcwd()
    candidates += [cwd, os.path.join(cwd, "F2_human_audit")]
    for directory in candidates:
        if os.path.exists(os.path.join(directory, sentinel)):
            return directory
    raise RuntimeError("cannot locate the audit .sage files; looked in: %s" % candidates)


AUDIT_DIR = _locate_audit_dir()
AUDIT_START_TIME = time.time()


def log(message):
    print("[audit %6.1fs] %s" % (time.time() - AUDIT_START_TIME, message))
    sys.stdout.flush()


def parse_command_line():
    parser = argparse.ArgumentParser(
        description="Audit the F2 genus-95 curve-code product property."
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="skip the exact symbolic audit and run only the random smoke test",
    )
    parser.add_argument(
        "--symbolic-only",
        action="store_true",
        help="run the exact symbolic audit and skip the random smoke test",
    )
    parser.add_argument(
        "--random-trials",
        type=int,
        default=RANDOM_TRIALS,
        help="number of random smoke-test trials to run",
    )
    # parse_known_args (not parse_args) so a launcher that leaves extra tokens on
    # the command line does not abort the audit.
    args, _ = parser.parse_known_args()

    # The new sage.cli launcher owns the command line and rejects script flags,
    # so also honor environment variables:
    #   AUDIT_QUICK=1, AUDIT_SYMBOLIC_ONLY=1, AUDIT_RANDOM_TRIALS=N
    if os.environ.get("AUDIT_QUICK", "") not in ("", "0"):
        args.quick = True
    if os.environ.get("AUDIT_SYMBOLIC_ONLY", "") not in ("", "0"):
        args.symbolic_only = True
    if os.environ.get("AUDIT_RANDOM_TRIALS"):
        args.random_trials = int(os.environ["AUDIT_RANDOM_TRIALS"])

    if args.quick and args.symbolic_only:
        parser.error("--quick and --symbolic-only cannot both be set")
    return args


if __name__ in ("__main__", "sage.all"):
    log("loading R and the product evaluator...")

load(os.path.join(AUDIT_DIR, "R_matrix_F2_genus95_candidate.sage"))
extension_matrix_R = R
product_coordinate_profile = list(row_profile)

load(os.path.join(AUDIT_DIR, "F2_genus95_candidate_product_evaluator_polynomials.sage"))
product_code_evaluate = eval_product_codeword

# Preserve the optional evaluator added upstream. Exact agreement is checked
# by certify.py; the legacy entry point also retains its runtime checks.
load(os.path.join(AUDIT_DIR, "F2_genus95_candidate_base_evaluator_polynomials.sage"))
base_code_evaluate_direct = eval_base_codeword

if product_coordinate_profile != list(product_row_profile):
    raise RuntimeError("R and the product evaluator use different coordinates")

if __name__ in ("__main__", "sage.all"):
    log("loaded %s product coordinates" % len(product_coordinate_profile))


# Explicit affine model of the curve over F2, in variables (x, y, z1, z2, z3).
# This is the single human-readable source for "what curve is this"; the
# coefficient lists below are exactly these equations and nothing else.
#
#   y^4 + (x+1)*y^2 + (x^3+x)*y + x^7 + x^3 = 0
#
#   z1^2 + z1 = [ (x^3+1)*y^3
#                 + (x^4+x^3+x^2)*y^2
#                 + (x^5+x^4+x^3+x^2+1)*y
#                 + (x^9+x^6+x^5+x^4+x^3+x^2) ] / (x^10+x^4+1)
#
#   z2^2 + z2 =   (x^4+x^2+x)/(x^11+x^10+x^5+x^4+x+1) * y^3
#               + [ (x^3+x)*y^2
#                   + (x^6+x^5+x^2+x)*y
#                   + (x^8+x^5+x^4+x^3+x^2+x) ] / (x^10+x^4+1)
#
#   z3^2 + z3 = [ (x^3+x^2)*y^3
#                 + (x^3+x^2+x+1)*y^2
#                 + (x^5+1)*y
#                 + (x^6+x^4+x^3+x^2) ] / (x^10+x^4+1)
#
# Note: (x^11+x^10+x^5+x^4+x+1 = (x+1)*(x^10+x^4+1), so z2's common denominator is
# (x+1)*(x^10+x^4+1).)
#
# Each equation is monic in a fresh variable -- degree 4 in y, degree 2 in each
# z_i -- so the curve is a degree-32 cover of the affine x-line: it is a curve.
# The coefficient lists are ordered [y^0, y^1, y^2, y^3].
artin_schreier_rhs_coeffs = [
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
artin_schreier_rhs = [
    _base_function_from_coeffs(coeffs)
    for coeffs in artin_schreier_rhs_coeffs
]
product_denominator_vector = [product_denominator] + [C(0)] * 7


def eval_rational_x_function(f, x0):
    K0 = parent(x0)
    num = f.numerator()
    den = f.denominator()
    den_value = K0(den(x0))
    if den_value == 0:
        raise ZeroDivisionError("base rational function denominator vanishes")
    return K0(num(x0)) / den_value


def eval_base_function(f, P):
    x0, y0 = P[0], P[1]
    K0 = parent(x0)
    total = K0(0)
    coeffs = list(f.list()) + [K(0)] * 4
    for i, coeff in enumerate(coeffs[:4]):
        if coeff:
            total += eval_rational_x_function(coeff, x0) * y0^i
    return total


row_index = {tuple(row): i for i, row in enumerate(product_coordinate_profile)}


def extend_base_message(base_message):
    message = vector(base_message) if hasattr(base_message, "base_ring") else vector(F2, base_message)
    if message.base_ring().characteristic() != 2:
        raise ValueError("messages must have characteristic two")
    return extension_matrix_R.change_ring(message.base_ring()) * message


def base_code_evaluate(base_message, P):
    """Single-source base evaluator: evaluate the product code on R*m."""
    return product_code_evaluate(extend_base_message(base_message), P)


def product_rule(left_product_message, right_product_message):
    """Multiply product-code coordinate vectors by the expected Hasse rules."""
    field = left_product_message.base_ring() if hasattr(left_product_message, "base_ring") else F2
    product_message = [field(0)] * len(product_coordinate_profile)

    for row in product_coordinate_profile:
        kind, point_index, jet_order = row
        out_index = row_index[tuple(row)]
        total = field(0)
        for left_order in range(jet_order + 1):
            right_order = jet_order - left_order
            left_index = row_index[(kind, point_index, left_order)]
            right_index = row_index[(kind, point_index, right_order)]
            total += (
                left_product_message[left_index]
                * right_product_message[right_index]
            )
        product_message[out_index] = total

    return vector(field, product_message)


# -----------------------------------------------------------------------------
# Symbolic reduced-function arithmetic
# -----------------------------------------------------------------------------

cover_basis = list(cover_basis_exponents)
cover_basis_index = {tuple(exps): i for i, exps in enumerate(cover_basis)}


def zero_reduced_function():
    return [C(0)] * len(cover_basis)


def reduced_monomial(exps, coefficient=1):
    out = zero_reduced_function()
    out[cover_basis_index[tuple(exps)]] = C(coefficient)
    return out


def add_reduced_functions(*functions):
    out = zero_reduced_function()
    for f in functions:
        for i, value in enumerate(f):
            out[i] += value
    return out


def reduce_monomial_product(exps1, exps2):
    terms = {
        tuple(exps1[i] + exps2[i] for i in range(3)): C(1)
    }

    changed = True
    while changed:
        changed = False
        new_terms = {}
        for exps, coefficient in terms.items():
            replacement_terms = None
            for variable_index in range(3):
                if exps[variable_index] >= 2:
                    changed = True
                    lowered_once = list(exps)
                    lowered_once[variable_index] -= 1
                    lowered_twice = list(exps)
                    lowered_twice[variable_index] -= 2
                    replacement_terms = [
                        (tuple(lowered_once), coefficient),
                        (
                            tuple(lowered_twice),
                            coefficient * artin_schreier_rhs[variable_index],
                        ),
                    ]
                    break

            if replacement_terms is None:
                replacement_terms = [(exps, coefficient)]

            for new_exps, new_coefficient in replacement_terms:
                new_terms[new_exps] = (
                    new_terms.get(new_exps, C(0)) + new_coefficient
                )

        terms = {
            exps: coefficient
            for exps, coefficient in new_terms.items()
            if coefficient
        }

    out = zero_reduced_function()
    for exps, coefficient in terms.items():
        out[cover_basis_index[exps]] += coefficient
    return out


monomial_products = None


def prepare_symbolic_arithmetic(show_progress=True):
    global monomial_products

    if monomial_products is not None:
        return

    if show_progress:
        log("preparing symbolic Artin-Schreier basis arithmetic...")

    monomial_products = {
        (i, j): reduce_monomial_product(cover_basis[i], cover_basis[j])
        for i in range(len(cover_basis))
        for j in range(len(cover_basis))
    }

    if show_progress:
        log("symbolic arithmetic table ready")


def multiply_reduced_functions(left, right):
    if monomial_products is None:
        prepare_symbolic_arithmetic(show_progress=False)

    out = zero_reduced_function()
    for i, left_coefficient in enumerate(left):
        if not left_coefficient:
            continue
        for j, right_coefficient in enumerate(right):
            if not right_coefficient:
                continue
            coefficient = left_coefficient * right_coefficient
            for k, product_coefficient in enumerate(monomial_products[(i, j)]):
                if product_coefficient:
                    out[k] += coefficient * product_coefficient
    return out


def combine_product_evaluator_numerators(product_message):
    out = zero_reduced_function()
    for coefficient, gamma in zip(product_message, gamma_numerators):
        if not coefficient:
            continue
        for i, gamma_coefficient in enumerate(gamma):
            if gamma_coefficient:
                out[i] += C(coefficient) * gamma_coefficient
    return out


def reduced_functions_are_equal(left, right):
    return all(left[i] == right[i] for i in range(len(cover_basis)))


# -----------------------------------------------------------------------------
# Nondegeneracy/rank checks
# -----------------------------------------------------------------------------
#
# These checks make the distance/injectivity assumptions executable.
#
#   - R must have rank 64, so distinct base messages give distinct product
#     coordinate vectors.
#   - The 222 exported numerator functions gamma_i must be F-linearly
#     independent over F = GF(2^128) -- that is what makes a nonzero F-message w
#     give a nonzero codeword Gamma_w = sum_i w_i gamma_i.  To check this
#     exactly, write every gamma_i in the 8 cover monomials and the 4 base
#     monomials 1,y,y^2,y^3, expand the F2(x) coefficients into their x-power
#     coefficients, and take the rank of the resulting matrix over F.  (The
#     gamma_i are F2-valued, so this equals the F2-rank; but we compute over F so
#     the check is the property actually needed, not one that happens to be
#     equivalent.)


def polynomial_function_coordinates(reduced_function):
    entries = []
    for base_coefficient in reduced_function:
        coefficients = list(base_coefficient.list()) + [K(0)] * 4
        entries.extend(coefficients[:4])
    return entries


def product_evaluator_numerator_matrix():
    # F-linear independence of the gamma_i (F = GF(2^128)) is what (A) needs, so
    # build the coefficient matrix over F.  Entries are the gamma_i F2(x)
    # coefficients expanded into x-power coefficients (each in F2 <= F).
    deploy_field = GF(2^EXTENSION_DEGREE_OVER_F2)
    polynomial_entries_by_gamma = [
        polynomial_function_coordinates(gamma)
        for gamma in gamma_numerators
    ]

    max_degree = 0
    for entries in polynomial_entries_by_gamma:
        for entry in entries:
            if entry.denominator() != 1:
                raise RuntimeError("product numerator entry is not polynomial")
            numerator = entry.numerator()
            if numerator:
                max_degree = max(max_degree, numerator.degree())

    columns = []
    for entries in polynomial_entries_by_gamma:
        flat_entries = []
        for entry in entries:
            numerator = entry.numerator()
            for degree in range(max_degree + 1):
                flat_entries.append(deploy_field(numerator[degree]))
        columns.append(vector(deploy_field, flat_entries))

    return matrix(deploy_field, columns).transpose()


def check_nondegeneracy(show_progress=True):
    if show_progress:
        log("checking rank(R)=64...")
    R_rank = extension_matrix_R.rank()
    if R_rank != 64:
        raise RuntimeError("R has rank %s, expected 64" % R_rank)

    if show_progress:
        log("checking product evaluator numerator rank=222...")
    evaluator_rank = product_evaluator_numerator_matrix().rank()
    expected_rank = len(product_coordinate_profile)
    if evaluator_rank != expected_rank:
        raise RuntimeError(
            "product evaluator numerator rank is %s, expected %s"
            % (evaluator_rank, expected_rank)
        )

    return R_rank, evaluator_rank


# -----------------------------------------------------------------------------
# Symbolic check
# -----------------------------------------------------------------------------
#
# The evaluator has the form
#
#     product_code_evaluate(w, P) = Gamma_w(P) / H(P),
#
# where Gamma_w is linear in the product-code coordinates w.  So the identity
#
#     evaluate(R*m, P) * evaluate(R*n, P)
#         = evaluate(product_rule(R*m, R*n), P)
#
# is equivalent, after multiplying by H(P)^2, to
#
#     Gamma_{R*m} * Gamma_{R*n}
#         = H * Gamma_{product_rule(R*m, R*n)}.
#
# The difference is bilinear in the two base messages.  Checking all pairs of
# standard basis vectors e_i,e_j therefore proves the identity for every pair
# of base messages.


def check_coordinate_reduction_rules():
    """Check that the reduction rules encode exactly the cover equations."""
    one = (0, 0, 0)
    generators = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

    for i, exps in enumerate(generators):
        expected = add_reduced_functions(
            reduced_monomial(exps),
            reduced_monomial(one, artin_schreier_rhs[i]),
        )
        actual = reduce_monomial_product(exps, exps)
        if not reduced_functions_are_equal(actual, expected):
            raise RuntimeError("coordinate reduction rule mismatch at z%s" % (i + 1))


def check_symbolic_basis_products(show_progress=True):
    if show_progress:
        log("checking coordinate reduction rules...")
    check_coordinate_reduction_rules()
    prepare_symbolic_arithmetic(show_progress=show_progress)

    if show_progress:
        log("building numerators for the 64 extended basis messages...")
    extended_basis_messages = [
        extension_matrix_R.column(i)
        for i in range(64)
    ]
    basis_numerators = [
        combine_product_evaluator_numerators(extended)
        for extended in extended_basis_messages
    ]

    checked = 0
    total = len(extended_basis_messages)^2
    if show_progress:
        log("checking %s symbolic basis-product identities..." % total)
    for i, left_extended in enumerate(extended_basis_messages):
        for j, right_extended in enumerate(extended_basis_messages):
            left_side = multiply_reduced_functions(
                basis_numerators[i],
                basis_numerators[j],
            )

            product_message = product_rule(left_extended, right_extended)
            product_numerator = combine_product_evaluator_numerators(product_message)
            right_side = multiply_reduced_functions(
                product_denominator_vector,
                product_numerator,
            )

            if not reduced_functions_are_equal(left_side, right_side):
                raise RuntimeError("symbolic product mismatch at (%s, %s)" % (i, j))
            checked += 1
            if (
                show_progress
                and (checked % SYMBOLIC_PROGRESS_INTERVAL == 0 or checked == total)
            ):
                log("symbolic basis products checked: %s/%s" % (checked, total))

    return checked


# -----------------------------------------------------------------------------
# Random check
# -----------------------------------------------------------------------------
#
# This is a smoke test for the same map-level identity at random affine curve
# points over GF(2^128).  It is not the proof; the symbolic check above is the
# exact audit.

RandomField.<b> = F2.extension(EXTENSION_DEGREE_OVER_F2)
RandomPoly.<T> = PolynomialRing(RandomField)


def random_base_message():
    return vector(F2, [F2.random_element() for _ in range(64)])


def random_product_message():
    return vector(
        F2,
        [F2.random_element() for _ in range(len(product_coordinate_profile))],
    )


def random_nonzero_product_message():
    for _ in range(MAX_POINT_ATTEMPTS):
        message = random_product_message()
        if any(coefficient != 0 for coefficient in message):
            return message
    raise RuntimeError("could not generate a nonzero product message")


def artin_schreier_roots(c):
    return (T^2 + T + RandomField(c)).roots(multiplicities=False)


def choose_root(roots):
    return roots[ZZ.random_element(len(roots))]


# -----------------------------------------------------------------------------
# What the sampler samples
# -----------------------------------------------------------------------------
#
# On each attempt choose x uniformly, reject d(x)*eta(x)=0, then choose
# one of FOUR slots for the ordered distinct y roots (reject missing slots).
# Each solvable Artin-Schreier equation has two roots, chosen uniformly.
# Thus every allowed point has proposal probability 1/(32*|F|); conditioning
# on acceptance is uniform on the paper's evaluation set Y.  The attempt cap
# raises an error rather than returning a biased fallback.


def random_curve_point():
    """Sample an affine point on the genus-95 cover over GF(2^128)."""
    for _ in range(MAX_POINT_ATTEMPTS):
        x0 = RandomField.random_element()
        if (x0 + 1) * (x0^10 + x0^4 + 1) == 0:
            continue
        if eval_base_polynomial_function(product_denominator, (x0, RandomField(0))) == 0:
            continue

        base_polynomial = (
            T^4
            + (x0 + 1) * T^2
            + (x0^3 + x0) * T
            + x0^7
            + x0^3
        )
        y_roots = base_polynomial.roots(multiplicities=False)
        y_roots.sort()
        slot = ZZ.random_element(4)
        if slot >= len(y_roots):
            continue
        y0 = y_roots[slot]

        z_values = []
        good_point = True
        for rhs in artin_schreier_rhs:
            rhs_value = eval_base_function(rhs, (x0, y0))
            roots = artin_schreier_roots(rhs_value)
            if not roots:
                good_point = False
                break
            z_values.append(choose_root(roots))
        if not good_point:
            continue

        return (x0, y0) + tuple(z_values)

    raise RuntimeError("could not find a random curve point")


def random_point_where_evaluator_is_defined():
    zero_product_message = vector(F2, [F2(0)] * len(product_coordinate_profile))

    for _ in range(MAX_POINT_ATTEMPTS):
        P = random_curve_point()
        try:
            product_code_evaluate(zero_product_message, P)
            return P
        except ZeroDivisionError:
            pass

    raise RuntimeError("could not find a point where the evaluator is defined")


def check_product_identity_at(P, left_base_message, right_base_message):
    left_extended = extend_base_message(left_base_message)
    right_extended = extend_base_message(right_base_message)

    process_1 = (
        base_code_evaluate(left_base_message, P)
        * base_code_evaluate(right_base_message, P)
    )

    product_message = product_rule(left_extended, right_extended)
    process_2 = product_code_evaluate(product_message, P)

    if process_1 != process_2:
        raise RuntimeError("product-property mismatch")


def run_random_checks(trials=RANDOM_TRIALS, show_progress=True):
    if show_progress:
        log("running %s random GF(2^128) smoke-test trials..." % trials)
    for trial in range(trials):
        P = random_point_where_evaluator_is_defined()
        check_product_identity_at(P, random_base_message(), random_base_message())
        if show_progress:
            log("random smoke-test trials checked: %s/%s" % (trial + 1, trials))
    return trials


def check_random_nonzero_codeword_evals(trials=RANDOM_TRIALS, show_progress=True):
    product_message = random_nonzero_product_message()
    if show_progress:
        log(
            "checking one random nonzero product codeword at %s random points..."
            % trials
        )

    for trial in range(trials):
        P = random_point_where_evaluator_is_defined()
        value = product_code_evaluate(product_message, P)
        if value == 0:
            raise RuntimeError(
                "random nonzero product codeword vanished at trial %s"
                % (trial + 1)
            )
        if show_progress:
            log(
                "nonzero codeword random evals checked: %s/%s"
                % (trial + 1, trials)
            )

    return trials


def check_base_evaluator_matches_product(trials=RANDOM_TRIALS, show_progress=True):
    """Cross-check the direct base evaluator against the product path.

    base_code_evaluate is defined as C(m)(P) := C_*(R*m)(P).  The direct
    evaluator computes the same value from the 64 base coordinates alone.  This
    checks they agree at random GF(2^128) points, i.e. that the standalone base
    evaluator file is a faithful shortcut.
    """
    if show_progress:
        log("cross-checking the direct base evaluator against C_*(R*m)...")
    for trial in range(trials):
        P = random_point_where_evaluator_is_defined()
        m = random_base_message()
        direct = base_code_evaluate_direct(m, P)
        via_product = base_code_evaluate(m, P)
        if direct != via_product:
            raise RuntimeError(
                "direct base evaluator disagrees with C_*(R*m) at trial %s"
                % (trial + 1)
            )
        if show_progress:
            log("base evaluator cross-checks: %s/%s" % (trial + 1, trials))
    return trials


if __name__ in ("__main__", "sage.all"):
    args = parse_command_line()
    log("starting audit")

    R_rank, evaluator_rank = check_nondegeneracy()
    print("R_rank:", R_rank)
    print("product_evaluator_numerator_rank:", evaluator_rank)
    print("nondegeneracy_check: ok")

    if args.quick:
        print("product_property_symbolic_check: skipped (--quick)")
    else:
        checked_symbolic_products = check_symbolic_basis_products()
        print("symbolic_basis_product_checks:", checked_symbolic_products)
        print("product_property_symbolic_check: ok")

    if args.symbolic_only:
        print("random_product_checks: skipped (--symbolic-only)")
        print("random_nonzero_codeword_evals: skipped (--symbolic-only)")
    else:
        print("field: GF(2^%s)" % EXTENSION_DEGREE_OVER_F2)
        checked_random_products = run_random_checks(trials=args.random_trials)
        print("random_product_checks:", checked_random_products)
        print("product_property_map_level_check: ok")
        checked_nonzero_evals = check_random_nonzero_codeword_evals(
            trials=args.random_trials
        )
        print("random_nonzero_codeword_evals:", checked_nonzero_evals)
        print("nonzero_codeword_random_eval_check: ok")
        checked_base_evaluator = check_base_evaluator_matches_product(trials=args.random_trials)
        print("base_evaluator_cross_checks:", checked_base_evaluator)
        print("base_evaluator_matches_product_check: ok")

    log("audit finished")


# -----------------------------------------------------------------------------
# Distance certificate
# -----------------------------------------------------------------------------
# Run python3 verify.py for the complete deterministic certificate, including
# the degree-32 field/full-constant-field check and the equation-only norm and
# point-count bounds.  See F2_human_audit/README.md for the proof and its mapping
# to the paper.  Running this product script alone does not certify distance.
