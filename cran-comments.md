## Test environments

* local: macOS, R 4.6.0
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1), macOS-latest
  (release), windows-latest (release)

## R CMD check results

0 errors | 0 warnings | 0 notes (aside from environment-specific notes such as
the inability to verify CRAN incoming feasibility when run offline).

## Notes for the maintainer

* This is a fork-maintained reliability release (0.3.0). The package bundles the
  archived `norm2` Fortran engine (Joseph L. Schafer / U.S. Census Bureau); its
  attribution is recorded in `inst/COPYRIGHTS`.
* `seed` now defaults to `NULL`, and the previously non-functional
  individual-specific-methods path (`methodvar`/`referencevar`) now errors
  explicitly rather than returning incomplete results. Both changes are described
  in `NEWS.md`.
