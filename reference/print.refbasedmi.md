# Print and summarise RefBasedMI results

[`RefBasedMI()`](https://choxos.github.io/RefbasedMI/reference/RefBasedMI.md)
returns its stacked long-format imputations as a data frame of class
`"refbasedmi"` carrying the run settings as attributes.
[`print()`](https://rdrr.io/r/base/print.html) shows a one-line
description of the run before the data;
[`summary()`](https://rdrr.io/r/base/summary.html) reports the settings
together with the row count and the number of missing outcome values in
each imputed dataset (which should be zero for every `.imp > 0`).

## Usage

``` r
# S3 method for class 'refbasedmi'
print(x, ...)

# S3 method for class 'refbasedmi'
summary(object, ...)
```

## Arguments

- x, object:

  An object returned by
  [`RefBasedMI()`](https://choxos.github.io/RefbasedMI/reference/RefBasedMI.md).

- ...:

  Passed on to the underlying data-frame method.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument
invisibly. [`summary()`](https://rdrr.io/r/base/summary.html) returns an
object of class `"summary.refbasedmi"`, a list with components `method`,
`reference`, `M`, `K0`, `K1`, `delta`, `dlag`, `depvar`, `n_rows`,
`rows_per_imputation`, and `missing_depvar_per_imputation`.

## Details

Subsetting keeps the class and the stored settings. If the attributes
have been removed, both methods fall back to the plain data-frame
behaviour.

## See also

[`RefBasedMI()`](https://choxos.github.io/RefbasedMI/reference/RefBasedMI.md),
[`as_mids()`](https://choxos.github.io/RefbasedMI/reference/as_mids.md)

## Examples

``` r
# \donttest{
asthmaJ2R <- RefBasedMI(data = asthma, depvar = fev, treatvar = treat,
                        idvar = id, timevar = time, covar = base,
                        method = "J2R", reference = 1, M = 2, seed = 101,
                        burnin = 100, verbose = FALSE)
summary(asthmaJ2R)
#> RefBasedMI imputed data
#>   Method:      J2R
#>   Reference:   1
#>   Imputations: M = 2
#>   Rows:        2196 (long format, .imp = 0 is the original data)
#> 
#> Missing 'fev' values by imputation number:
#>   0   1   2 
#> 147   0   0 
head(asthmaJ2R)
#> RefBasedMI imputed data: method J2R, reference 1, M = 2
#> Long format, 6 rows (.imp = 0 holds the original data)
#> 
#>      id time .imp  base treat  fev
#> 4  5001    2    0 2.925     1 2.87
#> 7  5001    4    0 2.925     1 2.66
#> 10 5001    8    0 2.925     1 2.69
#> 1  5001   12    0 2.925     1 2.58
#> 16 5003    2    0 2.465     2 2.61
#> 19 5003    4    0 2.465     2 2.45
# }
```
