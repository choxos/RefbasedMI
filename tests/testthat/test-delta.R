test_that("zero delta leaves imputations unchanged and non-zero delta shifts them", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  base_run <- quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                           timevar = time, covar = base, method = "J2R", reference = 1,
                           M = 2, seed = 5, burnin = 30)
  zero_delta <- quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                             timevar = time, covar = base, method = "J2R", reference = 1,
                             M = 2, seed = 5, burnin = 30,
                             delta = c(0, 0, 0, 0), dlag = c(1, 1, 1, 1))
  pos_delta <- quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                            timevar = time, covar = base, method = "J2R", reference = 1,
                            M = 2, seed = 5, burnin = 30,
                            delta = c(1, 1, 1, 1), dlag = c(1, 1, 1, 1))

  expect_equal(zero_delta$fev, base_run$fev)
  expect_false(isTRUE(all.equal(pos_delta$fev, base_run$fev)))
  # a positive delta can only raise post-discontinuation imputations
  expect_gte(sum(pos_delta$fev, na.rm = TRUE), sum(base_run$fev, na.rm = TRUE))
})
