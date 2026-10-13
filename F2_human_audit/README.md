# Certificate for the paper's explicit binary-preserving code

Read [**Speeding up Sumcheck using Algebraic Geometry** (PDF)](paper.pdf).

Run `python3 verify.py` from the repository root. The complete certificate
uses only the files in this directory, SageMath, and the root runner. It
establishes the equation-level claims in **Independent validation of the
explicit code**, without using the curve genus, a divisor, Riemann–Roch
spaces, or successful search logs as inputs.

## Claim-to-check map

The labels below are the paper's LaTeX labels, so they remain useful when
section and equation numbers change.

| Paper claim | Executable check | Expected result |
| --- | --- | --- |
| Printed matrix and `eq:eta` | `certify.py:check_paper_data` | All 222 rows and eta agree |
| `R = (I_64; R_0)`, block convolution | `certify.py:check_structure` | Binary 222-by-64 matrix; exact 64+64+64+30 order; formal polynomial identity for `B` |
| Product-systematicity and binary preservation | Same structure check | Zeroth jets multiply coordinatewise; all coefficients lie in GF(2) |
| `dim H = 222`, injective message embedding | `check_F2_product_properties.sage:check_nondegeneracy` | R rank 64; numerator rank 222 |
| `eq:degree-32-certificate`, full constant field | `check_F2_coordinate_field.sage:certify_coordinate_algebra_is_a_field`, called by `certify.py` | Degree 32; irreducible after extending constants to GF(2^128); full constants GF(2) |
| `eq:localized-equation-algebra` | Monic triangular presentation plus field check; denominator checks in `check_bounds` | Free monomial basis of size 32 over GF(2)[x,1/d] |
| Optional direct base evaluator | `certify.py:check_direct_base_evaluator` | All 64 rational-function identities agree with `C_*(R*v)` |
| `eq:productcheck` | `check_symbolic_basis_products` | All 4096 exact identities pass |
| Raw norm bound | `certify.py:check_bounds` | Cleared matrix weighted degree 80; determinant degree at most 2560; raw zero bound 83840 |
| `eq:explicit-norm-bound` | Same bounds check | Finite-pole clearing degree 16; norm numerator degree at most 2248 |
| `eq:explicit-point-count` | Same bounds check | Delta=13104; q>4 Delta^2; n >= 2^128-2^92 |
| `ssec:explicit-sampler` | `random_curve_point`; seeded runtime checks | Uniform bounded-fiber rejection, with d*eta nonzero |
| Independent relative distance | Integer inequalities in `check_bounds` | 2248/(2^128-2^92) < 2^-116 |

The certificate relies on exact Sage arithmetic and its polynomial
irreducibility and matrix-rank algorithms. The conclusions also use the
mathematical arguments below and the lemmas in the paper; this is not a
formal proof in a proof assistant. Six seeded evaluations exercise the
implementation but do not replace any symbolic check.

## Data and conventions

- `R_matrix_F2_genus95_candidate.sage` stores `R` and `row_profile`.
  Column `j` is bit `j` of a printed u64 row, with bit 0 least significant.
  The first 64 rows are the identity. The coordinate blocks are 64 values,
  64 first Hasse derivatives, 64 second Hasse derivatives, and the third
  Hasse derivatives at the first 30 points. All indices in the code are
  zero-based; `selected_third_derivative_indices_1_based` is explicitly
  one-based.
- `F2_genus95_candidate_product_evaluator_polynomials.sage` stores eta and
  the numerator table. It loads the same `row_profile` as R. The cover
  monomials are `(1,z3,z2,z2*z3,z1,z1*z3,z1*z2,z1*z2*z3)`.
  **Storage order is `[coordinate][cover monomial][y exponent]`**,
  whereas the paper writes `E[coordinate][y exponent][cover monomial]`.
  Thus `E[i,j,a] = gamma_numerator_coeffs[i][a][j]`.
- `F2_genus95_candidate_base_evaluator_polynomials.sage` provides the
  optional direct evaluator with 64 numerator functions and its own
  denominator. The certificate cross-multiplies denominators to check all
  64 basis functions against the defining product-evaluator path, and checks
  that its denominator has no roots on Y.
- `paper_instance.json` is a transcription of the printed 222 u64 rows and
  eta, with the paper's source revision. `--paper` additionally checks a
  supplied TeX source directly. The polynomial table is an input to the
  certificate, not a trusted result from a previous execution.

The denominator is

```text
eta = x^49+x^47+x^45+x^37+x^35+x^31+x^23+x^21+x^9+x^7+x^5+x.
```

The base equation is

```text
y^4 + (x+1)y^2 + (x^3+x)y + x^7+x^3 = 0.
```

The three equations `z_i^2+z_i=g_i` are printed as the coefficient lists
`artin_schreier_rhs_coeffs` in `check_F2_product_properties.sage`.
`certify.py` checks that the field certificate uses exactly the same lists
and base equation. Set `d=(x+1)(x^10+x^4+1)`. The evaluation set Y consists
of the solutions over GF(2^128) with `d(x)*eta(x) != 0`.
The auxiliary inverse coordinate in the paper is determined by x.

The maps are

```text
R(v) = R*v
C_*(w)(P) = sum_i w_i * gamma_i(P) / eta(P)
C(v)(P) = C_*(R*v)(P).
```

The implementation accepts Sage vectors over GF(2) or the deployment field.
Products are the truncated Hasse convolution independently at each point.
The zeroth block is therefore the pointwise product of the input messages.
The product check clears eta^2 and compares reduced functions for every pair
of the 64 base basis vectors. Bilinearity extends those equalities to all
messages over every extension field where evaluation is defined.

## Why the field check proves the point-count hypothesis

The monic tower has degrees 4,2,2,2, giving the 32 monomials
`y^j*z1^a*z2^b*z3^c`, with `0<=j<4` and `a,b,c` in `{0,1}`.
The field script forms multiplication by
`alpha=y+z1+x*z2+x^2*z3`, computes its degree-32 characteristic polynomial,
checks that it annihilates alpha, and tests irreducibility over GF(2)(x)
and GF(2^128)(x). Its irreducibility makes it the minimal polynomial and
identifies the explicit quotient with a field.

This also proves **geometric irreducibility**, rather than only
irreducibility over the deployment field. The full constant field has
degree c dividing 32. If c>1, then gcd(c,128)>1, and extending constants to
GF(2^128) would split the field. The deployment-field irreducibility test
rules this out. Thus c=1. The ground field is perfect, so this is the full
constant-field hypothesis in the paper. This inference specifically uses
the degrees 32 and 128; it is not true for arbitrary extensions.

A singular rational point alone does not establish this hypothesis.
The degree-32 statement is over GF(2)(x), not over GF(2).
After inverting d, successive division by the monic equations gives a
free rank-32 module. Its injection into its generic fiber follows from
freeness. This supplies the localized equation algebra and the bound of
32 points per fiber used in the norm argument.

## Recomputing the distance bounds

Write `p_w=sum_i w_i*gamma_i`. The rank-222 check makes this nonzero for
every nonzero w, including w over GF(2^128). The field check then makes its
norm nonzero. The norm-multiplicity lemma in the paper bounds the number
of zeros on Y by the numerator degree of that norm; the norm denominator
has roots only where d=0, which is excluded from Y.

For the raw bound, give x,y,z1,z2,z3 weights `1,7/4,0,0,0`.
The code checks the inequalities for the defining equations and inspects
every coefficient of E. In fact `degree_x E[i,j,a] <= 49-2*j`, so p_w has
weighted degree at most 49. The polynomial

```text
c_raw = (x^10+x^4+1)^3*(x+1)
```

of degree 31 clears the multiplication matrices. The certificate checks
all 32-by-32 basis reductions. If omega_i is the weight of basis monomial i,
the entry (i,j) has degree at most `80+omega_j-omega_i`. The weights cancel
in every determinant term, giving degree at most `32*80=2560`.
Individual entries can exceed degree 80: it is the **weighted** bound that
proves the determinant estimate. The paper's conservative fiber accounting
then gives `32*2560+32*(49+11)=83840` zeros.

For the sharper bound, put `h=x^5+x^2+1`, so `h^2=x^10+x^4+1`.
The script checks that

```text
h*z1, (x+1)*h*z2, h*z3
```

are integral over GF(2)[x,y]: after scaling, each satisfies a monic quadratic
with polynomial coefficients. Therefore `c*p_w` is integral over GF(2)[x],
where `c=(x+1)*h^3` has degree 16. The same assertions hold after extending
constants. The equations give pole weights at infinity at most `7/4` for
y and 0 for each z_i. The paper's conservative coefficient bound gives

```text
pole weight(c*p_w) <= 16 + 49 + 3*(7/4) = 281/4
norm numerator degree <= 32*(281/4) = 2248.
```

Here the norm of the integral element is a polynomial. Dividing it by
`c^32` cannot increase the reduced numerator degree. The checked exact
support permits a smaller estimate, but 2248 reproduces the stated bound.

For the evaluation-set size, clearing the four equations separately gives
total degrees 7,12,13,12, whose product is Delta=13104. The certificate
checks the hypotheses and numerical substitution in
[Cafure–Matera, Theorem 7.1](https://arxiv.org/html/math/0405302#S7), as used
in the paper. The absolutely irreducible component comes from the field
and full-constant-field checks. Removing d=0 costs at most 11*Delta points;
removing eta=0 afterwards costs at most 49*32 points. The script bounds
Delta^(13/3) upward using integer arithmetic, yielding

```text
n >= 2^128-2^92 > 2248
83840/n < 2^-111
2248/n < 2^-116.
```

No genus-95 point count is used. The first inequality also proves that
C_* is injective. The 316 zero bound and 119-bit claim require the
additional assertion that these functions represent sections of the
claimed degree-316 line bundle. The independent certificate does not
establish that assertion.

## Uniform sampler

`random_curve_point` first rejects roots of d*eta. It orders the distinct
y roots and chooses one of **four** slots uniformly, rejecting missing
slots. Each solvable Artin–Schreier equation then has exactly two roots;
the sampler chooses one uniformly at each of the three stages.
Every point in Y has proposal probability `1/(32*2^128)`, so the accepted
point is uniform. Its expected number of attempts is at most
`32*2^128/(2^128-2^92)`. The implementation imposes an attempt limit and
raises an error on exhaustion; it never returns a biased fallback.

## Tests of the checker

```sh
sage -python F2_human_audit/test_certificate.py
```

These tests deliberately corrupt a matrix entry, the paper denominator,
systematicity, the product prefix, the direct base evaluator, and evaluator independence, and require
the corresponding checks to reject. They also exercise the binary-list
interface. They test the checker; `python3 verify.py` runs the certificate.
