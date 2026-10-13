# AG-skip verification certificate

[![Verify certificate](https://github.com/succinctlabs/AG-skip-Verification-Certificate/actions/workflows/verify.yml/badge.svg)](https://github.com/succinctlabs/AG-skip-Verification-Certificate/actions/workflows/verify.yml)

An executable SageMath certificate for the explicit binary-preserving,
product-systematic code in **Speeding up Sumcheck using Algebraic Geometry**,
by Rahul Dalal, Benedikt Bünz, Tamir Hemo, and Ron D. Rothblum.

Read the [paper (PDF)](F2_human_audit/paper.pdf) and the
[detailed audit guide](F2_human_audit/README.md).

## Verify

Install [SageMath](https://www.sagemath.org/download.html) and Python 3.8 or
newer. The certificate is tested with **SageMath 10.8**; all mathematical
dependencies are included in SageMath. No separate Python packages are needed.

```sh
git clone https://github.com/succinctlabs/AG-skip-Verification-Certificate.git
cd AG-skip-Verification-Certificate
python3 verify.py
```

Allow about five minutes on a recent laptop; runtime varies by machine.
A successful full run exits with code 0 and ends with:

```text
result: PASS: equation-level certificate, 116-bit bound
```

The runner writes `.audit/verification.log` and `.audit/verification.json`.
The JSON report records the SageMath version, runtime, results, and SHA-256
hashes of the certificate inputs. These generated files are ignored by Git.
A failed check exits with a nonzero code; consult the log for details.

If `sage` is not on your PATH, use
`python3 verify.py --sage /path/to/sage`.
For a quick installation check, use `python3 verify.py --smoke`.
**Smoke mode skips the field and symbolic product certificates and does not
certify the result.**

## What is checked

The full run checks the printed matrix and denominator, the systematic
64-to-222 embedding and binary coefficients, the Hasse convolution product
rule, evaluator independence, all 4096 symbolic basis-product identities,
and all 64 identities for the direct base evaluator. It certifies that the
displayed degree-32 equation algebra is a field with full constant field
GF(2), and checks the bounds giving

```text
n >= 2^128 - 2^92
relative distance >= 1 - 2248/n > 1 - 2^-116.
```

This establishes the paper's independent **116-bit bound**. The **119-bit
geometric bound** also requires the curve and divisor interpretation and is
not established by this certificate. The audit guide maps each claim to its
check and explains the mathematical implications. The verification relies
on exact SageMath arithmetic and those arguments; it is not a proof-assistant
formalization. Seeded runtime checks also exercise the sampler and evaluators.

## Repository contents

| Path | Purpose |
| --- | --- |
| [`verify.py`](verify.py) | One-command runner and log capture |
| [`F2_human_audit/`](F2_human_audit/) | Self-contained certificate, input data, audit guide, tests, and paper |
| [`.github/workflows/verify.yml`](.github/workflows/verify.yml) | Full verification and rejection tests in a pinned SageMath 10.8 environment |

The certificate data and checking code come from revision
[`0e95326e83ab46a29828904bd5cadb1e7f916552`](https://github.com/succinctlabs/AG_codes/tree/0e95326e83ab46a29828904bd5cadb1e7f916552)
of the research repository. This release extracts the equation-field check
into `check_F2_coordinate_field.sage` and preserves the certificate data.
It includes no search runs, construction scripts, or F₄ material.

## Tests and paper cross-check

```sh
sage -python F2_human_audit/test_certificate.py
```

These seven tests check rejection of corrupted inputs and the binary-list
interface. They complement the full certificate. GitHub Actions runs both
commands and saves the verification report and log as a workflow artifact.

If you have the paper's TeX source, you can also compare its printed matrix
and denominator directly:

```sh
python3 verify.py --paper /path/to/explicit_instance.tex
```

Without this option, the certificate compares them against the included
`paper_instance.json` transcription. See [CONTRIBUTING.md](CONTRIBUTING.md)
for reporting problems or proposing changes.

## Citation

Rahul Dalal, Benedikt Bünz, Tamir Hemo, and Ron D. Rothblum.
*Speeding up Sumcheck using Algebraic Geometry*. 2026.
[Paper](F2_human_audit/paper.pdf) · [Citation metadata](CITATION.cff).
When reporting a verification run, also record the repository commit and
retain its JSON report.

## License

Dual-licensed under either of [Apache License, Version 2.0](LICENSE-APACHE)
or [MIT license](LICENSE-MIT) at your option, following
[Veil's formal verification repository](https://github.com/succinctlabs/veil-formal-verification#license).
