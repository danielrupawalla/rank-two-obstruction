# No four-fold unbiased measurements with four outcomes of rank two: verification package

Daniel Rupawalla, Independent Researcher (danielrupawalla@utexas.edu). 10 October 2026.

Repository: https://github.com/danielrupawalla/rank-two-obstruction

`rank_two_obstruction.tex` is the manuscript source. The scripts check the exact
algebraic identities used in its proofs. Case coverage and exhaustion are
analytic arguments in the manuscript; the scripts check the identities those
arguments rely on, not exhaustion by enumeration.

## Run

Requires Python 3.10+ and SymPy 1.14.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-verification.txt
.venv/bin/python run_manuscript_verifiers.py
```

The runner executes the fifteen verifiers listed in the manuscript's
verification appendix, with assertions enabled and a 180-second limit each.
Any failure or timeout gives a nonzero exit status. Each script writes a JSON
report; the runner writes `RUN_RESULTS.json`. Total runtime is under a minute
on a laptop.

## Files

| Manuscript result | Verifier |
|---|---|
| Moment identities (M4), (M6) and spectral heredity | `verify_spectral_moment_rigidity.py` |
| Triangle normal form | `verify_triangle_normal_form.py` |
| Scalar-triangle exclusion over F_3 | `verify_scalar_triangle_exclusion.py` |
| Commuting-pair nonextension | `verify_commuting_pair_theorem.py` |
| Fixed real pair: structure and nonextension | `verify_fixedpair_structure.py`, `verify_fixedpair_no_extension.py` |
| Coarse-support lemma and pair chart | `verify_coarse_support_bridge.py`, `verify_two_column_kernel.py`, `verify_unified_pair_chart.py` |
| Stratum symmetries | `verify_stratum_symmetry.py` |
| Global commuting-pair theorem | `verify_global_outcome_classification.py`, `verify_global_linear_reduction.py`, `verify_global_exceptional_reduction.py`, `verify_global_exceptional_complex_term.py`, `verify_global_coarse_reduction.py` |

Supporting files:

- `helpers.py`: shared exact quaternion-matrix arithmetic.
- `verify_extension_certificate.py`, `verify_family_reductions.py`, `verify_rotated_triple.py`: exact constructions imported by the fixed-pair and triangle verifiers.
- `extension_linear_certificate.json`: exact certificate read by `verify_extension_certificate.py`.
- `verify_mixed_strong_triple.py` (optional, not in the runner): an exact complete strong triple whose first pair is noncommuting while its other two pairs commute, illustrating the global commuting-pair theorem.

`MANIFEST.json` records SHA-256 hashes of the distributed files.

## License

Code: MIT (see `LICENSE`). The manuscript (`rank_two_obstruction.tex`, `.pdf`) is not covered by the code license; all rights reserved by the author.
