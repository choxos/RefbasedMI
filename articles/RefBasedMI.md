# 1. Getting started with RefBasedMI

## Background

In a longitudinal clinical trial, patients who deviate from the protocol
(for example by discontinuing randomised treatment) often have missing
outcomes afterwards. A primary analysis under “missing at random” (MAR)
assumes those patients would have continued along their own arm’s
trajectory. **Reference-based multiple imputation** instead encodes
assumptions in which deviating patients resemble a *reference* arm —
typically control — which is usually more plausible and more
conservative for an active treatment.

`RefBasedMI` implements these methods, following Carpenter, Roger and
Kenward (2013) and the causal model of White, Royes and Best (2020). It
is an R port of the Stata program *mimix*.

The reference-based methods differ in how the post-deviation mean is
constructed:

| Method | Assumption after deviation | Reference arm |
|----|----|----|
| MAR | continues on own arm | not required |
| J2R | jumps immediately to the reference mean | required |
| CR | follows the reference arm throughout | required |
| CIR | keeps the gained increment, then accrues reference increments | required |
| LMCF | last observed mean carried forward | not required |
| Causal | maintains a decaying fraction of the treatment effect | required |

## Quick start

The bundled `asthma` trial has the outcome `fev` measured at weeks 2, 4,
8 and 12, with a baseline covariate `base`. We impute under
jump-to-reference with control (arm 1) as reference, creating `M = 5`
completed datasets.

``` r

imp <- impute(
  data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
  covar = base, method = "J2R", reference = 1, M = 5, seed = 101
)
dim(imp)
#> [1] 4392    6
head(imp[, c(".imp", "id", "treat", "time", "fev")])
#>    .imp   id treat time  fev
#> 7     0 5001     1    2 2.87
#> 18    0 5001     1    4 2.66
#> 23    0 5001     1    8 2.69
#> 3     0 5001     1   12 2.58
#> 34    0 5003     2    2 2.61
#> 41    0 5003     2    4 2.45
```

The result is in long format: the original data (`.imp = 0`) stacked
above the `M` imputed datasets (`.imp = 1, ..., M`). Observed values are
unchanged; only the post-deviation missing values are filled in.
[`summary()`](https://rdrr.io/r/base/summary.html) reports the run
settings and confirms that every imputed dataset is complete:

``` r

summary(imp)
#> RefBasedMI imputed data
#>   Method:      J2R
#>   Reference:   1
#>   Imputations: M = 5
#>   Rows:        4392 (long format, .imp = 0 is the original data)
#> 
#> Missing 'fev' values by imputation number:
#>   0   1   2   3   4   5 
#> 147   0   0   0   0   0
```

### Analysis with Rubin’s rules

The output is designed to be turned into a `mids` object and pooled. The
exported
[`as_mids()`](https://choxos.github.io/RefbasedMI/reference/as_mids.md)
helper wraps
[`mice::as.mids()`](https://amices.org/mice/reference/as.mids.html):

``` r

fit <- with(as_mids(imp), lm(fev ~ factor(treat), subset = (time == 12)))
summary(mice::pool(fit))
#>             term  estimate std.error statistic       df      p.value
#> 1    (Intercept) 1.8864398 0.1048109 17.998506 17.77053 7.464593e-13
#> 2 factor(treat)2 0.2239927 0.1198322  1.869219 97.35187 6.459826e-02
```

## Comparing methods

Because the methods differ only in their assumption about deviating
patients, it is instructive to compare the imputed week-12 means across
methods using the same seed.

``` r

# RefBasedMI captures `method` by non-standard evaluation, so splice the literal
# method string into the call with bquote() rather than passing a loop variable.
impute_method <- function(m) {
  cl <- bquote(RefBasedMI(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = .(m), reference = 1, M = 5, seed = 101
  ))
  out <- NULL
  utils::capture.output(suppressMessages(out <- eval(cl)))
  out
}

methods <- c("MAR", "J2R", "CIR", "CR")
wk12 <- sapply(methods, function(m) {
  out <- impute_method(m)
  mean(out$fev[out$.imp > 0 & out$time == 12 & out$treat == 2])
})
round(wk12, 3)
#>   MAR   J2R   CIR    CR 
#> 2.208 2.110 2.166 2.182
```

The reference-based methods (J2R, CIR, CR) generally pull the active
arm’s imputed outcomes toward the control arm relative to MAR.

## The causal model

The causal model interpolates between J2R and CIR. The maintained
treatment effect after deviation is a fraction `K0` of the last
on-treatment effect, decaying geometrically by `K1` each period:
`K0 = 1, K1 = 0` reproduces J2R and `K0 = 1, K1 = 1` reproduces CIR.

``` r

acu <- impute(
  data = acupuncture, depvar = head, treatvar = treat, idvar = id, timevar = time,
  covar = head_base, method = "Causal", reference = 1, K0 = 1, K1 = 0.5,
  M = 5, seed = 101
)
dim(acu)
#> [1] 4812   12
```

## Delta adjustment

Delta adjustment adds a user-specified increment to post-deviation
imputations (but not to interim missing values), for tipping-point and
sensitivity analyses. `delta` gives the per-period increments and `dlag`
weights them. A constant shift of `-1` at the first post-deviation
visit, accruing thereafter, is:

``` r

impD <- impute(
  data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
  covar = base, method = "J2R", reference = 2,
  delta = c(-1, 0, 0, 0), dlag = c(1, 1, 1, 1), M = 5, seed = 101
)
dim(impD)
#> [1] 4392    6
```

## Interim (non-monotone) missingness

A missing value that is followed by a later observed value is treated as
*interim* and imputed under MAR regardless of the chosen reference-based
method; only post-discontinuation (monotone) missingness uses the
reference-based assumption. This matches the *mimix* behaviour and
requires no special handling from the user.

## Choosing M, burnin, and the seed

The examples in this vignette use `M = 5` (and `M = 2` in the help
pages) so they run quickly, but real analyses need more.

- **Number of imputations (`M`).** The Monte Carlo error of a pooled
  estimate falls as `M` grows. A common rule of thumb is to take `M` at
  least as large as the percentage fraction of missing information,
  which in trials with substantial post-deviation missingness often
  means 50 or more. Too small an `M` leaves the standard errors
  themselves noisy.
- **Burn-in (`burnin`) and spacing (`bbetween`).** Each arm’s
  multivariate normal parameters are drawn by the *norm2* MCMC sampler.
  `burnin` is the number of iterations discarded before the first draw;
  `bbetween` (defaulting to `burnin`) is the number of iterations
  between successive imputations, so that the `M` parameter draws are
  approximately independent. The defaults are generous for the small
  number of repeated measures typical of these trials; increase them if
  a trial has many time points or the sampler reports a near-boundary
  solution.
- **Reproducibility (`seed`).** `seed` defaults to `NULL`, so
  imputations are random each run. Pass an integer for a reproducible
  result.

A simple check that `M` is large enough is to repeat the whole
imputation under a few different seeds and confirm that the pooled
estimate barely moves:

``` r

wk12_estimate <- function(seed) {
  imp <- impute(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = "J2R", reference = 1, M = 20, seed = seed,
    verbose = FALSE
  )
  fit <- with(as_mids(imp), lm(fev ~ factor(treat), subset = (time == 12)))
  unname(mice::pool(fit)$pooled[["estimate"]][2])
}
round(vapply(c(101, 202, 303), wk12_estimate, numeric(1)), 3)
#> [1] 0.220 0.223 0.214
```

If those estimates differ materially, increase `M`.

## Other practical points

- **Baseline.** Supply the baseline outcome as a covariate (`covar`)
  rather than as an outcome value, so that no treatment effect is
  assumed at baseline. Factor covariates are accepted and expanded
  internally.
- **Silencing progress.** Pass `verbose = FALSE` to suppress the
  per-pattern progress messages, which is convenient inside loops such
  as the seed check above.
- **Individual-specific methods.** Assigning a different method per
  patient via `methodvar`/`referencevar` is not supported in this
  release. To vary the method across subgroups, impute each subgroup
  separately and combine the results.

## References

- Carpenter J R, Roger J H, Kenward M G (2013). Analysis of longitudinal
  trials with protocol deviation. *Journal of Biopharmaceutical
  Statistics* 23(6), 1352–1371.
- White I R, Royes J, Best N (2020). A causal modelling framework for
  reference-based imputation and tipping point analysis. *Journal of
  Biopharmaceutical Statistics* 30(2), 334–350.
