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
    out <- quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                        timevar = time, covar = base, method = cf$method,
                        reference = cf$reference, M = 2, seed = 1, burnin = 30)
    check_imputation_invariants(out, asthma, "fev", "id", "time", 2)
  }
})

test_that("the causal model yields a valid imputation", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  out <- quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                      timevar = time, covar = base, method = "Causal", reference = 1,
                      K0 = 1, K1 = 0.5, M = 2, seed = 1, burnin = 30)
  check_imputation_invariants(out, asthma, "fev", "id", "time", 2)
})
