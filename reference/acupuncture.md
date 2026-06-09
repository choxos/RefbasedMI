# Acupuncture trial data

Data from a randomised, double-blind, parallel-group trial comparing
active treatment with placebo. The outcome `head` (headache severity
score) is recorded in long format at times 3 and 12.

## Usage

``` r
acupuncture
```

## Format

A data frame with 802 rows and 11 variables:

- id:

  patient identifier

- time:

  measurement time point (3 or 12)

- age:

  patient age, in years

- sex:

  patient sex

- migraine:

  indicator of migraine (versus tension-type headache)

- chronicity:

  duration of headache disorder, in years

- practice_id:

  recruiting practice identifier

- treat:

  randomised treatment arm (1 or 2)

- head_base:

  baseline headache score, used as a covariate

- head:

  headache severity score, the outcome (missing after discontinuation)

- withdrawal_reason:

  reason for withdrawal, where recorded

## Source

Vickers A J et al. acupuncture for chronic headache trial.
