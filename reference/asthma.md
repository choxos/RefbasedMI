# Asthma trial data

The asthma clinical trial data used in the Stata *mimix* help file. The
outcome `fev` (forced expiratory volume) is recorded in long format at
weeks 2, 4, 8 and 12, with some post-baseline values missing.

## Usage

``` r
asthma
```

## Format

A data frame with 732 rows and 5 variables:

- id:

  patient identifier

- time:

  measurement time point, in weeks (2, 4, 8, 12)

- treat:

  randomised treatment arm (1 or 2)

- base:

  baseline value of the outcome, used as a covariate

- fev:

  forced expiratory volume, the outcome (missing after discontinuation)

## Source

Distributed with the Stata *mimix* package by Suzie Cro.
