# Reference-based multiple imputation of longitudinal clinical trial data

`RefBasedMI()` performs reference-based multiple imputation of a
continuous longitudinal outcome that is missing after treatment
discontinuation in a randomised clinical trial. It implements the
controlled multiple imputation methods of Carpenter, Roger and Kenward
(2013) – missing at random (MAR), jump to reference (J2R), copy
reference (CR), copy increments in reference (CIR) and last mean carried
forward (LMCF) – together with the causal model of White, Royes and Best
(2020) and delta adjustment for sensitivity analysis. It is an R port of
the Stata program *mimix* by Suzie Cro.

## Usage

``` r
RefBasedMI(
  data,
  covar = NULL,
  depvar,
  treatvar,
  idvar,
  timevar,
  method = NULL,
  reference = NULL,
  methodvar = NULL,
  referencevar = NULL,
  K0 = NULL,
  K1 = NULL,
  delta = NULL,
  dlag = NULL,
  M = 1,
  seed = NULL,
  prior = "jeffreys",
  burnin = 1000,
  bbetween = NULL,
  mle = FALSE,
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame of trial data in long format (one row per
  participant-time), sorted by participant.

- covar:

  Baseline covariate(s) given as a name or character vector; must be
  complete (no missing values) and numeric or factor. Typically the
  baseline outcome.

- depvar:

  Continuous outcome variable to be imputed.

- treatvar:

  Treatment arm variable; numeric or character.

- idvar:

  Participant identifier variable.

- timevar:

  Variable giving the repeated-measures time point.

- method:

  Group-level imputation method: one of `"MAR"`, `"J2R"`, `"CR"`,
  `"CIR"`, `"LMCF"` or `"Causal"` (case-insensitive). Exactly one of
  `method` or `methodvar` must be supplied.

- reference:

  Reference arm for `"J2R"`, `"CIR"`, `"CR"` and the control arm for
  `"Causal"`; numeric or character matching a level of `treatvar`.
  Required for those methods.

- methodvar:

  Variable specifying a per-participant method. Individual-specific
  methods are not supported in this release (see Details); supplying
  this raises an error.

- referencevar:

  Variable specifying a per-participant reference arm. Not supported in
  this release (see `methodvar`).

- K0:

  Causal model constant (the maintained fraction of the treatment
  effect); used with `method = "Causal"`.

- K1:

  Causal model decay constant in \[0, 1\]; the maintained effect decays
  by `K1` each period. Used with `method = "Causal"`.

- delta:

  Optional numeric vector of delta increments (the *a* values of Roger's
  "five macros"), of length equal to the number of time points.

- dlag:

  Optional numeric vector of delta lag weights (the *b* values), of
  length equal to the number of time points; defaults to all ones.

- M:

  Number of imputations to create.

- seed:

  Optional integer seed. Supply a value for reproducible imputations; if
  `NULL` (the default) the imputations are random each run.

- prior:

  Prior for the multivariate normal fit: `"jeffreys"` (default),
  `"uniform"` or `"ridge"`.

- burnin:

  Number of MCMC burn-in iterations.

- bbetween:

  Number of MCMC iterations between successive imputed datasets
  (defaults to `burnin`).

- mle:

  Logical; if `TRUE`, impute improperly by drawing from the maximum
  likelihood estimates instead of a posterior draw. Use with extreme
  caution: it ignores parameter uncertainty and invalidates
  Rubin's-rules intervals.

- verbose:

  Logical; if `TRUE` (the default) the function reports its progress,
  including the missing-data pattern summary and per-pattern imputation
  messages. Set to `FALSE` to run silently, which is convenient when
  calling `RefBasedMI()` repeatedly, for example over a sensitivity
  grid.

## Value

A data frame in long format stacking the original data (`.imp = 0`)
above the `M` imputed datasets (`.imp = 1, ..., M`). The `.imp` column
identifies the imputation and the participant identifier is retained.
Observed values are unchanged; only post-discontinuation (and interim)
missing outcomes are filled in. Pass the result to
[`mice::as.mids()`](https://amices.org/mice/reference/as.mids.html) to
analyse by Rubin's rules.

## Estimands and reference-based imputation

Patients who deviate from the protocol (most commonly by discontinuing
randomised treatment) often have unobserved post-deviation outcomes. The
analysis must state what those outcomes are assumed to be. A standard
MAR analysis assumes a deviating patient would have continued along
their own arm's trajectory. Reference-based methods instead assume the
patient comes to resemble a *reference* arm (usually control), which is
frequently more plausible and more conservative for an experimental
treatment.

## Algorithm

For each of `M` imputations the routine:

1.  builds a summary table of treatment arm by missing-data pattern;

2.  fits a multivariate normal distribution to each arm by MCMC (the
    bundled `norm2` engine);

3.  imputes interim (non-monotone) missing values under MAR;

4.  imputes post-discontinuation (monotone) missing values under the
    chosen `method`, constructing each patient's conditional mean from
    the relevant arm and reference-arm means;

5.  applies delta adjustment if requested.

The `M` completed datasets are stacked with the original data into one
long data frame for analysis by Rubin's rules.

## Methods

Writing the patient's own-arm and reference-arm mean trajectories as
\\\mu\\ and \\\mu^{ref}\\, the post-discontinuation conditional mean is:

- `"MAR"`:

  own-arm mean \\\mu\\; no reference required.

- `"J2R"`:

  jumps to the reference mean \\\mu^{ref}\\ from the first missing
  visit.

- `"CR"`:

  follows the reference mean \\\mu^{ref}\\ throughout.

- `"CIR"`:

  retains the treatment increment achieved up to deviation, then accrues
  reference-arm increments thereafter.

- `"LMCF"`:

  carries the last observed mean forward; no reference required.

- `"Causal"`:

  maintains a fraction of the last on-treatment effect that decays
  geometrically: the effect at lag \\k\\ is multiplied by \\K0 \times
  K1^{k}\\. `K0 = 1, K1 = 0` reproduces J2R and `K0 = 1, K1 = 1`
  reproduces CIR.

## Delta adjustment

Delta adjustment adds an increment to each post-discontinuation imputed
value (but not to interim missing values). For a patient who
discontinued after time \\p\\, the increment at time \\k\\ is
\\\delta\_{p+1} d_1 + \delta\_{p+2} d_2 + \dots + \delta_k d\_{k-p}\\,
where \\\delta\\ = `delta` and \\d\\ = `dlag`.

## Practical notes

Supply the baseline outcome as a covariate (`covar`) rather than as an
outcome value, so that no treatment effect is assumed at baseline.
Convert the output with
[`mice::as.mids()`](https://amices.org/mice/reference/as.mids.html) to
analyse it by Rubin's rules.

Individual-specific imputation methods (`methodvar` and `referencevar`)
are not supported in this release: the code path can leave records
unimputed and does not reproduce the corresponding group-level
imputations, so it raises an error rather than returning partial
results. To vary the method across subgroups, call `RefBasedMI()`
separately on each subgroup and combine the results.

## References

Carpenter J R, Roger J H, Kenward M G (2013). Analysis of longitudinal
trials with protocol deviation: a framework for relevant, accessible
assumptions, and inference via multiple imputation. *Journal of
Biopharmaceutical Statistics* 23(6), 1352–1371.
[doi:10.1080/10543406.2013.834911](https://doi.org/10.1080/10543406.2013.834911)

White I R, Royes J, Best N (2020). A causal modelling framework for
reference-based imputation and tipping point analysis in clinical trials
with quantitative outcome. *Journal of Biopharmaceutical Statistics*
30(2), 334–350.
[doi:10.1080/10543406.2019.1684308](https://doi.org/10.1080/10543406.2019.1684308)

## See also

[`mice::as.mids()`](https://amices.org/mice/reference/as.mids.html) and
[`mice::pool()`](https://amices.org/mice/reference/pool.html) for
analysing the output; the package vignettes
[`vignette("RefBasedMI")`](https://choxos.github.io/RefbasedMI/articles/RefBasedMI.md),
[`vignette("methods")`](https://choxos.github.io/RefbasedMI/articles/methods.md)
and
[`vignette("causal-and-delta")`](https://choxos.github.io/RefbasedMI/articles/causal-and-delta.md)
for worked examples.

## Examples

``` r
# Jump-to-reference imputation of the asthma trial, reference arm 1
asthmaJ2R <- RefBasedMI(data = asthma, depvar = fev, treatvar = treat,
  idvar = id, timevar = time, covar = base, method = "J2R", reference = 1,
  M = 2, seed = 54321)
#>    
#> 
#> Summary of missing data pattern by treat:
#> 
#>    pattern treat patients fev.2.miss fev.4.miss fev.8.miss fev.12.miss
#> 1        0     1       37          0          0          0           0
#> 2        7     1        1          1          1          1           0
#> 3        8     1       15          0          0          0           1
#> 4       12     1       22          0          0          1           1
#> 5       13     1        1          1          0          1           1
#> 6       14     1       16          0          1          1           1
#> 7        0     2       71          0          0          0           0
#> 8        4     2        1          0          0          1           0
#> 9        8     2        8          0          0          0           1
#> 10      12     2        8          0          0          1           1
#> 11      14     2        3          0          1          1           1
#> 
#> Fitting multivariate normal model by treat:
#>  
#> 
#> treat = 1
#> performing mcmcNorm for m = 1 to 2
#> 
#> mcmcNorm Loop finished.
#> 
#> treat = 2
#> performing mcmcNorm for m = 1 to 2
#> 
#> mcmcNorm Loop finished.
#> 
#> 
#> Number of original missing values = 147
#> 
#> Imputing interim missing values under MAR:
#> 
#> treat = 1pattern = 7number patients = 1
#> treat = 1pattern = 13number patients = 1
#> treat = 2pattern = 4number patients = 1
#> 
#> Number of post-discontinuation missing values = 142
#> 
#> Imputing post-discontinuation missing values under J2R:
#> 
#> treat  =  1 pattern =  0 number patients =  38 
#> treat  =  1 pattern =  8 number patients =  15 
#> treat  =  1 pattern =  12 number patients =  23 
#> treat  =  1 pattern =  14 number patients =  16 
#> treat  =  2 pattern =  0 number patients =  72 
#> treat  =  2 pattern =  8 number patients =  8 
#> treat  =  2 pattern =  12 number patients =  8 
#> treat  =  2 pattern =  14 number patients =  3 
#> 
#> Number of final missing values = 0
#> End of RefBasedMI

# Missing at random (no reference arm needed)
asthmaMAR <- RefBasedMI(data = asthma, depvar = fev, treatvar = treat,
  idvar = id, timevar = time, covar = base, method = "MAR", M = 2, seed = 54321)
#>    
#> 
#> Summary of missing data pattern by treat:
#> 
#>    pattern treat patients fev.2.miss fev.4.miss fev.8.miss fev.12.miss
#> 1        0     1       37          0          0          0           0
#> 2        7     1        1          1          1          1           0
#> 3        8     1       15          0          0          0           1
#> 4       12     1       22          0          0          1           1
#> 5       13     1        1          1          0          1           1
#> 6       14     1       16          0          1          1           1
#> 7        0     2       71          0          0          0           0
#> 8        4     2        1          0          0          1           0
#> 9        8     2        8          0          0          0           1
#> 10      12     2        8          0          0          1           1
#> 11      14     2        3          0          1          1           1
#> 
#> Fitting multivariate normal model by treat:
#>  
#> 
#> treat = 1
#> performing mcmcNorm for m = 1 to 2
#> 
#> mcmcNorm Loop finished.
#> 
#> treat = 2
#> performing mcmcNorm for m = 1 to 2
#> 
#> mcmcNorm Loop finished.
#> 
#> 
#> Number of original missing values = 147
#> 
#> Imputing interim missing values under MAR:
#> 
#> treat = 1pattern = 7number patients = 1
#> treat = 1pattern = 13number patients = 1
#> treat = 2pattern = 4number patients = 1
#> 
#> Number of post-discontinuation missing values = 142
#> 
#> Imputing post-discontinuation missing values under MAR:
#> 
#> treat  =  1 pattern =  0 number patients =  38 
#> treat  =  1 pattern =  8 number patients =  15 
#> treat  =  1 pattern =  12 number patients =  23 
#> treat  =  1 pattern =  14 number patients =  16 
#> treat  =  2 pattern =  0 number patients =  72 
#> treat  =  2 pattern =  8 number patients =  8 
#> treat  =  2 pattern =  12 number patients =  8 
#> treat  =  2 pattern =  14 number patients =  3 
#> 
#> Number of final missing values = 0
#> End of RefBasedMI

# Copy increments in reference
asthmaCIR <- RefBasedMI(data = asthma, depvar = fev, treatvar = treat,
  idvar = id, timevar = time, covar = base, method = "CIR", reference = 1,
  M = 2, seed = 54321)
#>    
#> 
#> Summary of missing data pattern by treat:
#> 
#>    pattern treat patients fev.2.miss fev.4.miss fev.8.miss fev.12.miss
#> 1        0     1       37          0          0          0           0
#> 2        7     1        1          1          1          1           0
#> 3        8     1       15          0          0          0           1
#> 4       12     1       22          0          0          1           1
#> 5       13     1        1          1          0          1           1
#> 6       14     1       16          0          1          1           1
#> 7        0     2       71          0          0          0           0
#> 8        4     2        1          0          0          1           0
#> 9        8     2        8          0          0          0           1
#> 10      12     2        8          0          0          1           1
#> 11      14     2        3          0          1          1           1
#> 
#> Fitting multivariate normal model by treat:
#>  
#> 
#> treat = 1
#> performing mcmcNorm for m = 1 to 2
#> 
#> mcmcNorm Loop finished.
#> 
#> treat = 2
#> performing mcmcNorm for m = 1 to 2
#> 
#> mcmcNorm Loop finished.
#> 
#> 
#> Number of original missing values = 147
#> 
#> Imputing interim missing values under MAR:
#> 
#> treat = 1pattern = 7number patients = 1
#> treat = 1pattern = 13number patients = 1
#> treat = 2pattern = 4number patients = 1
#> 
#> Number of post-discontinuation missing values = 142
#> 
#> Imputing post-discontinuation missing values under CIR:
#> 
#> treat  =  1 pattern =  0 number patients =  38 
#> treat  =  1 pattern =  8 number patients =  15 
#> treat  =  1 pattern =  12 number patients =  23 
#> treat  =  1 pattern =  14 number patients =  16 
#> treat  =  2 pattern =  0 number patients =  72 
#> treat  =  2 pattern =  8 number patients =  8 
#> treat  =  2 pattern =  12 number patients =  8 
#> treat  =  2 pattern =  14 number patients =  3 
#> 
#> Number of final missing values = 0
#> End of RefBasedMI

# Causal model: half the treatment effect retained, decaying by 0.5 per period
asthmaCausal <- RefBasedMI(data = asthma, depvar = fev, treatvar = treat,
  idvar = id, timevar = time, covar = base, method = "Causal", reference = 1,
  K0 = 1, K1 = 0.5, M = 2, seed = 54321)
#>    
#> 
#> Summary of missing data pattern by treat:
#> 
#>    pattern treat patients fev.2.miss fev.4.miss fev.8.miss fev.12.miss
#> 1        0     1       37          0          0          0           0
#> 2        7     1        1          1          1          1           0
#> 3        8     1       15          0          0          0           1
#> 4       12     1       22          0          0          1           1
#> 5       13     1        1          1          0          1           1
#> 6       14     1       16          0          1          1           1
#> 7        0     2       71          0          0          0           0
#> 8        4     2        1          0          0          1           0
#> 9        8     2        8          0          0          0           1
#> 10      12     2        8          0          0          1           1
#> 11      14     2        3          0          1          1           1
#> 
#> Fitting multivariate normal model by treat:
#>  
#> 
#> treat = 1
#> performing mcmcNorm for m = 1 to 2
#> 
#> mcmcNorm Loop finished.
#> 
#> treat = 2
#> performing mcmcNorm for m = 1 to 2
#> 
#> mcmcNorm Loop finished.
#> 
#> 
#> Number of original missing values = 147
#> 
#> Imputing interim missing values under MAR:
#> 
#> treat = 1pattern = 7number patients = 1
#> treat = 1pattern = 13number patients = 1
#> treat = 2pattern = 4number patients = 1
#> 
#> Number of post-discontinuation missing values = 142
#> 
#> Imputing post-discontinuation missing values under CAUSAL:
#> 
#> treat  =  1 pattern =  0 number patients =  38 
#> treat  =  1 pattern =  8 number patients =  15 
#> treat  =  1 pattern =  12 number patients =  23 
#> treat  =  1 pattern =  14 number patients =  16 
#> treat  =  2 pattern =  0 number patients =  72 
#> treat  =  2 pattern =  8 number patients =  8 
#> treat  =  2 pattern =  12 number patients =  8 
#> treat  =  2 pattern =  14 number patients =  3 
#> 
#> Number of final missing values = 0
#> End of RefBasedMI

# Delta adjustment: shift post-discontinuation imputations down by 1 unit
asthmaDelta <- RefBasedMI(data = asthma, depvar = fev, treatvar = treat,
  idvar = id, timevar = time, covar = base, method = "J2R", reference = 2,
  delta = c(-1, 0, 0, 0), dlag = c(1, 1, 1, 1), M = 2, seed = 54321)
#>    
#> 
#> Summary of missing data pattern by treat:
#> 
#>    pattern treat patients fev.2.miss fev.4.miss fev.8.miss fev.12.miss
#> 1        0     1       37          0          0          0           0
#> 2        7     1        1          1          1          1           0
#> 3        8     1       15          0          0          0           1
#> 4       12     1       22          0          0          1           1
#> 5       13     1        1          1          0          1           1
#> 6       14     1       16          0          1          1           1
#> 7        0     2       71          0          0          0           0
#> 8        4     2        1          0          0          1           0
#> 9        8     2        8          0          0          0           1
#> 10      12     2        8          0          0          1           1
#> 11      14     2        3          0          1          1           1
#> 
#> Fitting multivariate normal model by treat:
#>  
#> 
#> treat = 1
#> performing mcmcNorm for m = 1 to 2
#> 
#> mcmcNorm Loop finished.
#> 
#> treat = 2
#> performing mcmcNorm for m = 1 to 2
#> 
#> mcmcNorm Loop finished.
#> 
#> 
#> Number of original missing values = 147
#> 
#> Imputing interim missing values under MAR:
#> 
#> treat = 1pattern = 7number patients = 1
#> treat = 1pattern = 13number patients = 1
#> treat = 2pattern = 4number patients = 1
#> 
#> Number of post-discontinuation missing values = 142
#> 
#> Imputing post-discontinuation missing values under J2R:
#> 
#> treat  =  1 pattern =  0 number patients =  38 
#> treat  =  1 pattern =  8 number patients =  15 
#> treat  =  1 pattern =  12 number patients =  23 
#> treat  =  1 pattern =  14 number patients =  16 
#> treat  =  2 pattern =  0 number patients =  72 
#> treat  =  2 pattern =  8 number patients =  8 
#> treat  =  2 pattern =  12 number patients =  8 
#> treat  =  2 pattern =  14 number patients =  3 
#> 
#> Number of final missing values = 0
#> End of RefBasedMI

# Analyse the imputations with Rubin's rules via the mice package
# \donttest{
if (requireNamespace("mice", quietly = TRUE)) {
  fit <- with(mice::as.mids(asthmaJ2R),
              lm(fev ~ factor(treat), subset = (time == 12)))
  summary(mice::pool(fit))
}
#>             term  estimate  std.error statistic       df      p.value
#> 1    (Intercept) 1.8484905 0.09751412  18.95613 11.77294 3.480994e-10
#> 2 factor(treat)2 0.2586256 0.12345074   2.09497 66.75392 3.997138e-02
# }
```
