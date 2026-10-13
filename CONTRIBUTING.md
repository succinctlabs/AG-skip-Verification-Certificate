# Contributing

Please report questions, mismatches with the paper, and verification failures
in [GitHub Issues](https://github.com/succinctlabs/AG-skip-Verification-Certificate/issues).
For a failure, include the repository commit (`git rev-parse HEAD`), operating
system, SageMath version (`sage --version`), command, and relevant log output.
Attach `.audit/verification.json` when one was produced.

Before proposing a change, run:

```sh
sage -python F2_human_audit/test_certificate.py
python3 verify.py
```

Explain any effect on the certified claims and update the audit guide when
needed. Changes to certificate data must be justified against the paper;
update the paper transcription only when the printed instance changes.
Include a targeted rejection test when correcting a missed invalid input.
Keep generated reports, exploratory code, and construction outputs out of
this verification repository.

Contributions are provided under the repository's dual MIT/Apache-2.0 license.
