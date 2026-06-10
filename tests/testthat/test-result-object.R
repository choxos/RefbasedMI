impute_j2r_small <- function() {
  data(asthma, package = "RefBasedMI")
  quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
               timevar = time, covar = base, method = "J2R", reference = 1,
               M = 2, seed = 1, burnin = 30)
}

test_that("the result carries the refbasedmi class and run settings", {
  out <- impute_j2r_small()
  expect_s3_class(out, "refbasedmi")
  expect_s3_class(out, "data.frame")
  expect_identical(attr(out, "method"), "J2R")
  expect_identical(attr(out, "reference"), "1")
  expect_identical(attr(out, "M"), 2)
  expect_null(attr(out, "K0"))
  expect_null(attr(out, "delta"))
  expect_identical(attr(out, "depvar"), "fev")
  expect_identical(attr(out, "idvar"), "id")
})

test_that("print and summary report the run settings", {
  out <- impute_j2r_small()

  # subsetting keeps the class and settings, so the header survives head()
  printed <- capture.output(print(head(out, 2)))
  expect_true(any(grepl("method J2R, reference 1, M = 2", printed)))

  # a stripped copy falls back to plain data-frame printing
  bare <- out
  attr(bare, "method") <- NULL
  printed_bare <- capture.output(print(head(bare, 2)))
  expect_false(any(grepl("RefBasedMI imputed data", printed_bare)))

  s <- summary(out)
  expect_s3_class(s, "summary.refbasedmi")
  expect_identical(s$method, "J2R")
  expect_identical(as.vector(s$missing_depvar_per_imputation[c("1", "2")]),
                   c(0L, 0L))
  expect_true(s$missing_depvar_per_imputation[["0"]] > 0)
  shown <- capture.output(print(s))
  expect_true(any(grepl("Method:      J2R", shown)))
})

test_that("the result still behaves as a data frame", {
  out <- impute_j2r_small()
  expect_true(is.data.frame(out))
  sub <- out[out$.imp == 1, ]
  expect_s3_class(sub, "data.frame")
  expect_equal(nrow(sub), nrow(out) / 3)
})

test_that("as_mids round-trips into a mids object", {
  skip_if_not_installed("mice")
  out <- impute_j2r_small()
  mids <- as_mids(out)
  expect_s3_class(mids, "mids")
  expect_equal(mids$m, 2)
  pooled <- mice::pool(with(mids, lm(fev ~ factor(treat), subset = (time == 12))))
  expect_true(all(is.finite(summary(pooled)$estimate)))
})

test_that("summary attributes reflect causal and delta settings", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  out <- quiet_impute_method(asthma, method = "Causal", reference = 1,
                             K0 = 1, K1 = 0.5)
  expect_identical(attr(out, "method"), "Causal")
  expect_identical(attr(out, "K0"), 1)
  expect_identical(attr(out, "K1"), 0.5)

  outd <- quiet_impute(data = asthma, depvar = fev, treatvar = treat,
                       idvar = id, timevar = time, covar = base,
                       method = "J2R", reference = 1, M = 2, seed = 1,
                       burnin = 30, delta = c(0, 0.5, 1, 1))
  expect_identical(attr(outd, "delta"), c(0, 0.5, 1, 1))
  shown <- capture.output(print(summary(outd)))
  expect_true(any(grepl("Delta:       0, 0.5, 1, 1", shown)))
})
