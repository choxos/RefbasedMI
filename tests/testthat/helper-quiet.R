# RefBasedMI is intentionally chatty (progress messages mirroring the Stata mimix
# port). Wrap calls so the test output stays readable.
quiet_impute <- function(...) {
  out <- NULL
  suppressWarnings(suppressMessages(
    utils::capture.output(out <- RefBasedMI(...))
  ))
  out
}

# RefBasedMI captures `method` and `reference` by non-standard evaluation, so a
# value held in a variable must be spliced into the call as a literal (passing
# the variable directly would deparse the variable expression, not its value).
quiet_impute_method <- function(dat, method, reference = NULL, K0 = NULL, K1 = NULL,
                                M = 2, seed = 1, burnin = 30) {
  env <- new.env()
  env$dat <- dat
  cl <- bquote(RefBasedMI(
    data = dat, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = .(method), reference = .(reference),
    K0 = .(K0), K1 = .(K1), M = .(M), seed = .(seed), burnin = .(burnin)
  ))
  out <- NULL
  suppressWarnings(suppressMessages(
    utils::capture.output(out <- eval(cl, env))
  ))
  out
}
