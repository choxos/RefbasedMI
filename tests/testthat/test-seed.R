test_that("a supplied seed makes imputation reproducible", {
  data(asthma, package = "RefBasedMI")
  a <- quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                    timevar = time, covar = base, method = "J2R", reference = 1,
                    M = 2, seed = 123, burnin = 30)
  b <- quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                    timevar = time, covar = base, method = "J2R", reference = 1,
                    M = 2, seed = 123, burnin = 30)
  expect_equal(a, b)
})

test_that("different seeds give different imputations", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  a <- quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                    timevar = time, covar = base, method = "J2R", reference = 1,
                    M = 2, seed = 1, burnin = 30)
  b <- quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                    timevar = time, covar = base, method = "J2R", reference = 1,
                    M = 2, seed = 2, burnin = 30)
  expect_false(isTRUE(all.equal(a[a$.imp > 0, "fev"], b[b$.imp > 0, "fev"])))
})
