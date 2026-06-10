# Edge cases drawn from realistic trial data shapes.

test_that("character participant ids reproduce the numeric-id imputation", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  asthma_pid <- asthma
  # same sort order as the numeric ids, so the imputation must be identical
  asthma_pid$id <- paste0("P", asthma$id)

  out_chr <- quiet_impute(data = asthma_pid, depvar = fev, treatvar = treat,
                          idvar = id, timevar = time, covar = base,
                          method = "J2R", reference = 1, M = 2,
                          seed = 1, burnin = 30)
  out_num <- quiet_impute(data = asthma, depvar = fev, treatvar = treat,
                          idvar = id, timevar = time, covar = base,
                          method = "J2R", reference = 1, M = 2,
                          seed = 1, burnin = 30)
  expect_identical(out_chr$fev, out_num$fev)
  expect_identical(sort(unique(out_chr$id)), sort(unique(asthma_pid$id)))
  expect_false(any(is.na(out_chr$fev[out_chr$.imp > 0])))
})

test_that("a factor covariate matches its manual dummy expansion", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  fc <- asthma
  fc$site <- factor(ifelse(fc$id %% 2 == 0, "A", "B"))
  out_fac <- quiet_impute(data = fc, depvar = fev, treatvar = treat,
                          idvar = id, timevar = time, covar = c(base, site),
                          method = "J2R", reference = 1, M = 2,
                          seed = 1, burnin = 30)

  md <- fc
  md$site_B <- as.numeric(md$site == "B")
  out_man <- quiet_impute(data = md, depvar = fev, treatvar = treat,
                          idvar = id, timevar = time, covar = c(base, site_B),
                          method = "J2R", reference = 1, M = 2,
                          seed = 1, burnin = 30)

  expect_identical(out_fac$fev, out_man$fev)
  # the internal dummy column is not leaked into the output
  expect_false("site_B" %in% names(out_fac))
  expect_s3_class(out_fac$site, "factor")
})

test_that("a participant with every visit missing is fully imputed", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  allmiss <- asthma
  first_id <- allmiss$id[1]
  allmiss$fev[allmiss$id == first_id] <- NA

  out <- quiet_impute(data = allmiss, depvar = fev, treatvar = treat,
                      idvar = id, timevar = time, covar = base,
                      method = "J2R", reference = 1, M = 2, seed = 1,
                      burnin = 30)
  expect_false(any(is.na(out$fev[out$.imp > 0])))
  expect_equal(sum(out$id == first_id & out$.imp > 0), 4 * 2)
})

test_that("a single-arm trial can be imputed under MAR", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  one_arm <- asthma[asthma$treat == 1, ]
  out <- quiet_impute(data = one_arm, depvar = fev, treatvar = treat,
                      idvar = id, timevar = time, covar = base,
                      method = "MAR", M = 1, seed = 1, burnin = 30)
  expect_false(any(is.na(out$fev[out$.imp > 0])))
  expect_equal(nrow(out), nrow(one_arm) * 2)
})

test_that("missing participant-visit rows are completed and imputed", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  holes <- asthma[-c(2, 7, 20), ]
  out <- quiet_impute(data = holes, depvar = fev, treatvar = treat,
                      idvar = id, timevar = time, covar = base,
                      method = "J2R", reference = 1, M = 1, seed = 1,
                      burnin = 30)
  # the wide reshape restores the full participant-by-visit grid
  expect_equal(nrow(out), length(unique(holes$id)) * 4 * 2)
  expect_false(any(is.na(out$fev[out$.imp > 0])))
})

test_that("an all-zero delta reproduces the unadjusted imputation", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  d0 <- quiet_impute(data = asthma, depvar = fev, treatvar = treat,
                     idvar = id, timevar = time, covar = base,
                     method = "J2R", reference = 1, M = 1, seed = 3,
                     burnin = 30, delta = c(0, 0, 0, 0))
  pl <- quiet_impute(data = asthma, depvar = fev, treatvar = treat,
                     idvar = id, timevar = time, covar = base,
                     method = "J2R", reference = 1, M = 1, seed = 3,
                     burnin = 30)
  expect_identical(d0$fev, pl$fev)
})

test_that("time values starting at zero are handled", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  t0 <- asthma
  t0$time <- t0$time - 2
  out <- quiet_impute(data = t0, depvar = fev, treatvar = treat,
                      idvar = id, timevar = time, covar = base,
                      method = "J2R", reference = 1, M = 1, seed = 1,
                      burnin = 30)
  expect_setequal(unique(out$time), c(0, 2, 6, 10))
  expect_false(any(is.na(out$fev[out$.imp > 0])))
})
