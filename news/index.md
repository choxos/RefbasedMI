# Changelog

## RefBasedMI 0.3.0

This is a maintenance and reliability release. The reference-based
imputation results for the supported group-level methods (MAR, J2R, CR,
CIR, LMCF, Causal) are unchanged: every change below was verified
against a 14-scenario regression baseline captured from version 0.2.0,
and all scenarios remain numerically identical unless noted.

### Behaviour changes

- `seed` now defaults to `NULL` instead of `101`. Imputations are random
  by default and reproducible only when a seed is supplied. Code that
  relied on the implicit default must now pass `seed = 101` explicitly.
- Individual-specific imputation via `methodvar` / `referencevar` now
  raises an explicit error. The feature was non-functional in previous
  versions (it could leave records unimputed and did not reproduce the
  group-level imputations), so it is disabled rather than allowed to
  return partial results. To vary the method across subgroups, impute
  each subgroup separately with the group-level `method` and `reference`
  arguments and combine the results.

### Bug fixes

- Fixed native routine registration: the C initialisation routine was
  misnamed (a copy-paste leftover) and registered no routines, leaving
  the package reliant on dynamic symbol lookup. The five Fortran
  routines of the bundled `norm2` engine are now properly registered.
- `loglikNorm()` dispatched to the wrong S3 generic and computed the
  log-posterior instead of the log-likelihood; corrected.
- Treatment- and id-column handling used substring matching (`grep`) and
  boundary-unsafe index arithmetic; replaced with exact matching and a
  position-safe permutation.
- Corrected the delta guard `length(delta != 0)` and replaced
  `1:length(x)` loop bounds with
  [`seq_along()`](https://rdrr.io/r/base/seq.html) where the vector can
  be empty.
- Repaired several latent errors in the individual-methods entry path (a
  short-circuit logical bug, a missing non-standard-evaluation of
  `referencevar`, and a data-frame indexing error in the causal helper)
  ahead of the gate above.

### Dependencies

- Reduced runtime dependencies to base R (`stats`, `utils`). Removed
  `Hmisc` (replaced by an internal, order-identical grouping helper),
  `data.table`, `pastecs`, and `assertthat`. Moved `mice` to `Suggests`
  (used only in examples and vignettes).

### Documentation and infrastructure

- Rewrote the reference manual: a structured
  [`RefBasedMI()`](https://choxos.github.io/RefbasedMI/reference/RefBasedMI.md)
  help page with per-method descriptions, the delta formula, complete
  parameter documentation, references, and worked examples; complete
  dataset documentation.
- Added three vignettes: a getting-started guide, a methods comparison
  with imputed-trajectory plots, and a sensitivity-analysis vignette
  covering the causal model and delta tipping-point analysis.
- Added a `testthat` suite, GitHub Actions for `R CMD check` and
  `lintr`, and a pkgdown website.
- Added Ahmad Sofi-Mahmudi as an author.
