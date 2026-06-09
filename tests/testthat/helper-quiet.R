# RefBasedMI is intentionally chatty (progress messages mirroring the Stata mimix
# port). Wrap calls so the test output stays readable.
quiet_impute <- function(...) {
  out <- NULL
  suppressWarnings(suppressMessages(
    utils::capture.output(out <- RefBasedMI(...))
  ))
  out
}
