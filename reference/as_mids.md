# Convert a RefBasedMI result to a mice::mids object

A convenience wrapper around
[`mice::as.mids()`](https://amices.org/mice/reference/as.mids.html) so
the imputations can be analysed with
[`mice::with.mids()`](https://amices.org/mice/reference/with.mids.html)
and pooled by Rubin's rules with
[`mice::pool()`](https://amices.org/mice/reference/pool.html). Requires
the suggested package mice.

## Usage

``` r
as_mids(x)
```

## Arguments

- x:

  An object returned by
  [`RefBasedMI()`](https://choxos.github.io/RefbasedMI/reference/RefBasedMI.md)
  (or any stacked long-format data frame with an `.imp` column, where
  `.imp = 0` is the original data).

## Value

A [`mice::mids()`](https://amices.org/mice/reference/mids.html) object
holding the `M` imputed datasets.

## See also

[`RefBasedMI()`](https://choxos.github.io/RefbasedMI/reference/RefBasedMI.md),
[`mice::as.mids()`](https://amices.org/mice/reference/as.mids.html),
[`mice::pool()`](https://amices.org/mice/reference/pool.html)

## Examples

``` r
# \donttest{
if (requireNamespace("mice", quietly = TRUE)) {
  asthmaJ2R <- RefBasedMI(data = asthma, depvar = fev, treatvar = treat,
                          idvar = id, timevar = time, covar = base,
                          method = "J2R", reference = 1, M = 5, seed = 101,
                          burnin = 100, verbose = FALSE)
  fit <- with(as_mids(asthmaJ2R),
              lm(fev ~ factor(treat), subset = (time == 12)))
  summary(mice::pool(fit))
}
#>             term  estimate  std.error statistic       df      p.value
#> 1    (Intercept) 1.8682957 0.09217855 20.268227 40.88325 6.117558e-23
#> 2 factor(treat)2 0.2463132 0.12354486  1.993715 70.86655 5.003283e-02
# }
```
