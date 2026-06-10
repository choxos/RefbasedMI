test_that("missing or misspecified variables are rejected", {
  data(asthma, package = "RefBasedMI")

  expect_error(
    quiet_impute(data = asthma, depvar = notacolumn, treatvar = treat,
                 idvar = id, timevar = time, method = "MAR", M = 1),
    "depvar"
  )
  expect_error(
    quiet_impute(data = asthma, depvar = fev, treatvar = treat,
                 idvar = id, timevar = time, method = "NOTAMETHOD",
                 reference = 1, M = 1)
  )
  expect_error(
    quiet_impute(data = asthma, depvar = fev, treatvar = treat,
                 idvar = id, timevar = time, method = "J2R",
                 reference = 99, M = 1),
    "reference"
  )
})

test_that("reference-based methods require a reference arm", {
  data(asthma, package = "RefBasedMI")
  expect_error(
    quiet_impute(data = asthma, depvar = fev, treatvar = treat,
                 idvar = id, timevar = time, method = "J2R", M = 1),
    "reference"
  )
})

test_that("method and methodvar are mutually exclusive", {
  data(asthma, package = "RefBasedMI")
  expect_error(
    quiet_impute(data = asthma, depvar = fev, treatvar = treat,
                 idvar = id, timevar = time, method = "J2R",
                 methodvar = treat, reference = 1, M = 1),
    "NOT both"
  )
})

test_that("individual-specific methods are gated with a clear message", {
  data(asthma, package = "RefBasedMI")
  expect_error(
    quiet_impute(data = asthma, depvar = fev, treatvar = treat,
                 idvar = id, timevar = time,
                 methodvar = treat, referencevar = treat, M = 1),
    "not\\s+supported"
  )
})

test_that("scalar control arguments are validated early", {
  data(asthma, package = "RefBasedMI")
  j2r <- function(...) {
    quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                 timevar = time, method = "J2R", reference = 1, ...)
  }
  # a fractional M would silently truncate the parameter-draw index
  expect_error(j2r(M = 0), "M must be")
  expect_error(j2r(M = 2.5), "M must be")
  expect_error(j2r(M = c(2, 3)), "M must be")
  expect_error(j2r(M = 1, burnin = 0), "burnin must be")
  expect_error(j2r(M = 1, burnin = 10.5), "burnin must be")
  expect_error(j2r(M = 1, bbetween = -1), "bbetween must be")
  expect_error(j2r(M = 1, delta = c("a", "b", "c", "d")), "delta must be")
  expect_error(j2r(M = 1, delta = c(1, 1, 1, 1), dlag = letters[1:4]),
               "dlag must be")
  expect_error(j2r(M = 1, mle = c(TRUE, TRUE)), "mle must be")
})

test_that("data content problems are rejected before any model fitting", {
  data(asthma, package = "RefBasedMI")

  asthma_na <- asthma
  asthma_na$treat[1] <- NA
  expect_error(
    quiet_impute(data = asthma_na, depvar = fev, treatvar = treat, idvar = id,
                 timevar = time, method = "J2R", reference = 1, M = 1),
    "missing values in treat"
  )

  # a duplicated id-time row would be dropped silently by reshape()
  asthma_dup <- asthma[c(1, seq_len(nrow(asthma))), ]
  expect_error(
    quiet_impute(data = asthma_dup, depvar = fev, treatvar = treat, idvar = id,
                 timevar = time, method = "J2R", reference = 1, M = 1),
    "duplicate"
  )

  asthma_chr <- asthma
  asthma_chr$fev <- as.character(asthma_chr$fev)
  expect_error(
    quiet_impute(data = asthma_chr, depvar = fev, treatvar = treat, idvar = id,
                 timevar = time, method = "J2R", reference = 1, M = 1),
    "fev \\(depvar\\) must be numeric"
  )

  asthma_t <- asthma
  asthma_t$time <- as.character(asthma_t$time)
  expect_error(
    quiet_impute(data = asthma_t, depvar = fev, treatvar = treat, idvar = id,
                 timevar = time, method = "J2R", reference = 1, M = 1),
    "time \\(timevar\\) must be numeric"
  )
})

test_that("Causal constants are validated as numeric scalars", {
  data(asthma, package = "RefBasedMI")
  causal <- function(...) {
    quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                 timevar = time, method = "Causal", reference = 1, M = 1, ...)
  }
  expect_error(causal(), "K0 Causal constant not specified")
  # an explicit NULL used to crash with "argument is of length zero"
  expect_error(causal(K0 = NULL, K1 = 1), "K0 Causal constant not specified")
  expect_error(causal(K0 = "a", K1 = 1), "K0 Causal constant must be")
  expect_error(causal(K0 = 1), "K1 Causal constant not specified")
  expect_error(causal(K0 = 1, K1 = 2), "not in range")
})

test_that("K0/K1 supplied with a non-Causal method draw a warning", {
  data(asthma, package = "RefBasedMI")
  expect_warning(
    utils::capture.output(suppressMessages(
      RefBasedMI(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                 timevar = time, method = "J2R", reference = 1, K0 = 1,
                 M = 1, seed = 1, burnin = 10)
    )),
    "only used by the Causal method"
  )
})
