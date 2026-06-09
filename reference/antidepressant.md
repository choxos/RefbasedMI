# Antidepressant trial data

The antidepressant clinical trial data described by White, Royes and
Best (2020). The outcome `HAMD17.TOTAL` (Hamilton depression rating) is
recorded in long format at visits 4, 5, 6 and 7. The data also carry
per-patient method and reference columns used to illustrate
individual-specific imputation.

## Usage

``` r
antidepressant
```

## Format

A data frame with 688 rows and 14 variables:

- PATIENT.NUMBER:

  patient identifier

- HAMA.TOTAL:

  Hamilton anxiety rating scale total

- PGI.IMPROVEMENT:

  patient global impression of improvement

- VISIT.DATE:

  day of the visit relative to randomisation

- VISIT.NUMBER:

  visit number (4, 5, 6, 7)

- TREATMENT.NAME:

  randomised treatment arm

- PATIENT.SEX:

  patient sex

- POOLED.INVESTIGATOR:

  pooled investigator centre

- basval:

  baseline value of the outcome, used as a covariate

- HAMD17.TOTAL:

  Hamilton 17-item depression rating, the outcome

- change:

  change from baseline in the outcome

- miss_flag:

  indicator that the record is observed

- methodcol:

  per-patient imputation method (for individual-specific imputation)

- referencecol:

  per-patient reference arm (for individual-specific imputation)

## Source

White I R, Royes J, Best N (2020)
[doi:10.1080/10543406.2019.1684308](https://doi.org/10.1080/10543406.2019.1684308)
.
