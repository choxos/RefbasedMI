## Test environments

* local: macOS, R 4.6.0
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1), macOS-latest
  (release), windows-latest (release)

## R CMD check results

0 errors | 0 warnings | 0 notes (aside from environment-specific notes such as
the inability to verify CRAN incoming feasibility when run offline).

## Notes for the maintainer

* This is a fork-maintained release (0.3.1). The package bundles the archived
  `norm2` Fortran engine (Joseph L. Schafer / U.S. Census Bureau); its
  attribution is recorded in `inst/COPYRIGHTS`.
* 0.3.1 adds input validation, a `verbose` argument, a classed result object
  with `print()`/`summary()` methods, and an `as_mids()` helper. The
  reference-based imputation results are unchanged: every change was verified
  against a 14-scenario regression baseline and remains numerically identical.
  All changes are described in `NEWS.md`.
