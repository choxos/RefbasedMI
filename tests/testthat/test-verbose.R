test_that("verbose = FALSE silences all progress output", {
  data(asthma, package = "RefBasedMI")
  out <- NULL
  msgs <- character()
  stdout_lines <- withCallingHandlers(
    utils::capture.output(
      out <- RefBasedMI(data = asthma, depvar = fev, treatvar = treat,
                        idvar = id, timevar = time, covar = base,
                        method = "J2R", reference = 1, M = 1, seed = 1,
                        burnin = 30, verbose = FALSE)
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_identical(stdout_lines, character(0))
  expect_identical(msgs, character(0))
  expect_s3_class(out, "data.frame")
})

test_that("verbose = TRUE keeps the progress reporting", {
  data(asthma, package = "RefBasedMI")
  msgs <- testthat::capture_messages(
    utils::capture.output(
      RefBasedMI(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                 timevar = time, covar = base, method = "J2R", reference = 1,
                 M = 1, seed = 1, burnin = 30)
    )
  )
  expect_true(any(grepl("Fitting multivariate normal model", msgs)))
  expect_true(any(grepl("End of RefBasedMI", msgs)))
})

test_that("verbose must be a single logical", {
  data(asthma, package = "RefBasedMI")
  expect_error(
    quiet_impute(data = asthma, depvar = fev, treatvar = treat, idvar = id,
                 timevar = time, method = "J2R", reference = 1, M = 1,
                 verbose = "yes"),
    "verbose"
  )
})
