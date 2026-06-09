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
