# Character treatment-arm labels (e.g. "Placebo") exercise a different path
# from numeric-string labels: the original labels cannot be recovered by
# as.numeric(), so the internal 1..k codes are kept and the labels are
# reapplied at the end.

strip_settings <- function(d) {
  d <- as.data.frame(d)
  attributes(d) <- attributes(d)[c("names", "row.names", "class")]
  class(d) <- "data.frame"
  d
}

test_that("character arm labels reproduce the numeric-coded imputation", {
  data(asthma, package = "RefBasedMI")
  asthma_chr <- asthma
  # alphabetical order matches the numeric coding, so the same seed must give
  # exactly the same imputed values
  asthma_chr$treat <- c("active", "placebo")[asthma$treat]

  out_chr <- quiet_impute(data = asthma_chr, depvar = fev, treatvar = treat,
                          idvar = id, timevar = time, covar = base,
                          method = "J2R", reference = "placebo", M = 2,
                          seed = 1, burnin = 30)
  out_num <- quiet_impute(data = asthma, depvar = fev, treatvar = treat,
                          idvar = id, timevar = time, covar = base,
                          method = "J2R", reference = 2, M = 2,
                          seed = 1, burnin = 30)

  expect_identical(strip_settings(out_chr)$fev, strip_settings(out_num)$fev)
  expect_identical(attr(out_chr, "reference"), "placebo")

  # per-row arm labels must match the input data
  orig <- asthma_chr[order(asthma_chr$id, asthma_chr$time), ]
  imp0 <- out_chr[out_chr$.imp == 0, ]
  imp0 <- imp0[order(imp0$id, imp0$time), ]
  expect_identical(as.character(imp0$treat), as.character(orig$treat))
  expect_false(any(is.na(out_chr$fev[out_chr$.imp > 0])))
})

test_that("arm labels that flip the alphabetical order keep the mapping", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  flipped <- asthma
  flipped$treat <- c("zzz", "aaa")[asthma$treat]

  out <- quiet_impute(data = flipped, depvar = fev, treatvar = treat,
                      idvar = id, timevar = time, covar = base,
                      method = "J2R", reference = "aaa", M = 2,
                      seed = 1, burnin = 30)
  expect_false(any(is.na(out$fev[out$.imp > 0])))

  orig <- flipped[order(flipped$id, flipped$time), ]
  imp0 <- out[out$.imp == 0, ]
  imp0 <- imp0[order(imp0$id, imp0$time), ]
  expect_identical(as.character(imp0$treat), as.character(orig$treat))
})

test_that("a factor treatment variable with character levels works", {
  skip_on_cran()
  data(asthma, package = "RefBasedMI")
  fac <- asthma
  fac$treat <- factor(c("Drug", "Placebo")[asthma$treat])

  out <- quiet_impute(data = fac, depvar = fev, treatvar = treat,
                      idvar = id, timevar = time, covar = base,
                      method = "J2R", reference = "Placebo", M = 1,
                      seed = 1, burnin = 30)
  expect_s3_class(out$treat, "factor")
  expect_setequal(levels(out$treat), c("Drug", "Placebo"))
  expect_false(any(is.na(out$fev[out$.imp > 0])))
})
