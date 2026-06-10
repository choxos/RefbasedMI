# RefBasedMI 0.3.2

This release makes character treatment-arm labels, character participant ids,
and factor covariates work end to end, and adds guidance to the error messages.
The imputation results for previously supported inputs are unchanged (all 14
regression baseline scenarios remain byte-identical).

## Bug fixes

* Character treatment-arm labels (for example `"Placebo"` and drug names)
  previously crashed with `invalid 'labels'; length k should be 1 or 0` after
  `NAs introduced by coercion` warnings, and the missing-data pattern summary
  displayed the treatment column as `NA`. The internal numeric round-trip now
  runs only for arms whose labels parse as numbers; character labels keep their
  internal codes and are restored on output, and the pattern summary shows the
  real arm names.
* Character participant ids crashed in the interim-imputation machinery (data
  was round-tripped through `as.matrix()`, coercing every column to character).
  Ids are now recoded to integers internally and restored on output; with ids
  of equal sort order the imputed values are identical to a numeric-id run.
* Factor covariates, documented as supported, crashed in the conditional-draw
  arithmetic. They are now expanded into numeric dummy columns on entry
  (matching the `model.matrix()` design) and the helper columns are removed
  from the returned dataset; a factor covariate reproduces its manual dummy
  encoding exactly.
* `covar = c(...)` detection no longer depends on `sQuote()` fancy quotes,
  which are disabled in some environments (for example inside testthat).

## Improved error messages

* An invalid `reference` now lists the available treatment values; an unsorted
  id column shows the sorting one-liner; an unrecognised `method` lists the
  valid methods; covariate errors name the offending columns; `delta`/`dlag`
  length errors state the expected length; and the "covariance matrix is not
  positive definite" failure suggests `prior = "ridge"`, fewer covariates, or a
  longer burn-in.

## Tests

* New edge-case suites cover character arms (order-preserving and
  order-flipping labels, factor treatment variables), character ids, factor
  covariates, an all-visits-missing participant, single-arm MAR runs, missing
  participant-visit rows, all-zero delta equivalence, and time values starting
  at zero.

# RefBasedMI 0.3.1

This release adds input validation, a `verbose` option, and a richer result
object on top of 0.3.0. The reference-based imputation results for the supported
group-level methods are unchanged: every change was verified against the
14-scenario regression baseline and remains numerically identical.

## New features

* `RefBasedMI()` gains a `verbose` argument (default `TRUE`). Set `verbose =
  FALSE` to run silently, which is convenient when calling the function in a
  loop, for example over a delta sensitivity grid.
* The result is now a data frame of class `"refbasedmi"` carrying the run
  settings as attributes, with `print()` and `summary()` methods. `summary()`
  reports the settings together with the number of remaining missing outcome
  values in each imputed dataset.
* New exported `as_mids()` helper wraps `mice::as.mids()` so the imputations can
  be pooled by Rubin's rules in one call.

## Bug fixes

* The Causal model no longer errors when `K0` or `K1` is passed as an explicit
  `NULL`. When `K0 = 0` the decay constant `K1` may be omitted, as documented;
  this configuration now reproduces jump-to-reference exactly (previously it
  failed with "argument is of length zero").
* A non-integer `M` is rejected instead of silently truncating the
  parameter-draw index, which could select another arm's posterior draw in
  trials with three or more arms.
* The internally generated Fortran seeds are floored at 1, since `runif()` can in
  principle return 0.
* The `mle` argument is validated as a scalar, avoiding a length-mismatch error
  under R >= 4.2.

## Input validation

* `M`, `burnin`, and `bbetween` must be positive integers; `delta` and `dlag`
  must be numeric; `depvar` and `timevar` columns must be numeric.
* Missing values in the treatment, identifier, or time variables, and duplicated
  identifier-by-time rows (which `reshape()` would silently drop), are rejected
  at entry with a clear message.
* Supplying `K0` or `K1` with a non-Causal method now warns rather than being
  silently ignored.

## Documentation

* The getting-started vignette gains a section on choosing the number of
  imputations, the MCMC burn-in and spacing, and the seed, with a worked
  seed-sensitivity check.

# RefBasedMI 0.3.0

This is a maintenance and reliability release. The reference-based imputation
results for the supported group-level methods (MAR, J2R, CR, CIR, LMCF, Causal)
are unchanged: every change below was verified against a 14-scenario regression
baseline captured from version 0.2.0, and all scenarios remain numerically
identical unless noted.

## Behaviour changes

* `seed` now defaults to `NULL` instead of `101`. Imputations are random by
  default and reproducible only when a seed is supplied. Code that relied on the
  implicit default must now pass `seed = 101` explicitly.
* Individual-specific imputation via `methodvar` / `referencevar` now raises an
  explicit error. The feature was non-functional in previous versions (it could
  leave records unimputed and did not reproduce the group-level imputations), so
  it is disabled rather than allowed to return partial results. To vary the method
  across subgroups, impute each subgroup separately with the group-level `method`
  and `reference` arguments and combine the results.

## Bug fixes

* Fixed native routine registration: the C initialisation routine was misnamed
  (a copy-paste leftover) and registered no routines, leaving the package reliant
  on dynamic symbol lookup. The five Fortran routines of the bundled `norm2`
  engine are now properly registered.
* `loglikNorm()` dispatched to the wrong S3 generic and computed the log-posterior
  instead of the log-likelihood; corrected.
* Treatment- and id-column handling used substring matching (`grep`) and
  boundary-unsafe index arithmetic; replaced with exact matching and a
  position-safe permutation.
* Corrected the delta guard `length(delta != 0)` and replaced `1:length(x)` loop
  bounds with `seq_along()` where the vector can be empty.
* Repaired several latent errors in the individual-methods entry path (a
  short-circuit logical bug, a missing non-standard-evaluation of `referencevar`,
  and a data-frame indexing error in the causal helper) ahead of the gate above.

## Dependencies

* Reduced runtime dependencies to base R (`stats`, `utils`). Removed `Hmisc`
  (replaced by an internal, order-identical grouping helper), `data.table`,
  `pastecs`, and `assertthat`. Moved `mice` to `Suggests` (used only in examples
  and vignettes).

## Documentation and infrastructure

* Rewrote the reference manual: a structured `RefBasedMI()` help page with
  per-method descriptions, the delta formula, complete parameter documentation,
  references, and worked examples; complete dataset documentation.
* Added three vignettes: a getting-started guide, a methods comparison with
  imputed-trajectory plots, and a sensitivity-analysis vignette covering the
  causal model and delta tipping-point analysis.
* Added a `testthat` suite, GitHub Actions for `R CMD check` and `lintr`, and a
  pkgdown website.
* Added Ahmad Sofi-Mahmudi as an author.
