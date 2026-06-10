# Shared invariant checks for an imputed long-format result.
check_imputation_invariants <- function(out, data, depvar, idvar, timevar, M) {
  expect_s3_class(out, "data.frame")
  expect_true(".imp" %in% names(out))
  expect_setequal(unique(out$.imp), 0:M)
  expect_equal(nrow(out), nrow(data) * (M + 1))

  blocks <- split(out, out$.imp)
  ord <- function(b) b[order(b[[idvar]], b[[timevar]]), ]
  orig <- ord(blocks[["0"]])

  for (m in seq_len(M)) {
    imp <- ord(blocks[[as.character(m)]])
    # imputed sets are complete for the outcome
    expect_false(any(is.na(imp[[depvar]])),
                 info = paste("NA remains in imputation", m))
    # the same records are present
    expect_equal(imp[[idvar]], orig[[idvar]])
    expect_equal(imp[[timevar]], orig[[timevar]])
    # observed outcome values are left untouched
    obs <- !is.na(orig[[depvar]])
    expect_equal(imp[[depvar]][obs], orig[[depvar]][obs],
                 info = paste("observed values altered in imputation", m))
  }
}

test_that("J2R produces a complete, structure-preserving imputation (fast)", {
  data(asthma, package = "RefBasedMI")
  out <- quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                      timevar = time, covar = base, method = "J2R", reference = 1,
                      M = 2, seed = 1, burnin = 30)
  check_imputation_invariants(out, asthma, "fev", "id", "time", 2)
})

test_that("every reference-based method yields a valid imputation", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  cfg <- list(
    list(method = "MAR",  reference = NULL),
    list(method = "J2R",  reference = 1),
    list(method = "CR",   reference = 1),
    list(method = "CIR",  reference = 1),
    list(method = "LMCF", reference = NULL)
  )
  for (cf in cfg) {
    out <- quiet_impute_method(asthma, method = cf$method, reference = cf$reference)
    check_imputation_invariants(out, asthma, "fev", "id", "time", 2)
  }
})

test_that("the causal model yields a valid imputation", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  out <- quiet_impute_method(asthma, method = "Causal", reference = 1,
                             K0 = 1, K1 = 0.5)
  check_imputation_invariants(out, asthma, "fev", "id", "time", 2)
})

test_that("the causal model with K0 = 0 collapses to jump-to-reference", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  # K0 = 0 makes the maintained-effect term vanish, so K1 may be omitted and
  # the imputation must match J2R exactly under the same seed
  out0 <- quiet_impute_method(asthma, method = "Causal", reference = 1, K0 = 0)
  check_imputation_invariants(out0, asthma, "fev", "id", "time", 2)
  outj <- quiet_impute_method(asthma, method = "J2R", reference = 1)
  # the stored run settings differ by design; the imputed data must not
  strip_settings <- function(d) {
    d <- as.data.frame(d)
    attributes(d) <- attributes(d)[c("names", "row.names", "class")]
    class(d) <- "data.frame"
    d
  }
  expect_equal(strip_settings(out0), strip_settings(outj))
})
