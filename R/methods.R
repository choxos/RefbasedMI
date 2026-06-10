# S3 class, methods, and mice bridge for the RefBasedMI result object.

# Attach the refbasedmi class and run metadata to the assembled long-format
# imputation data. Called once, at the end of a successful RefBasedMI() run.
as_refbasedmi <- function(x, method, reference, M, K0, K1, delta, dlag,
                          depvar, treatvar, idvar, timevar) {
  attr(x, "method") <- method
  attr(x, "reference") <- reference
  attr(x, "M") <- M
  attr(x, "K0") <- if (length(K0) == 1) K0 else NULL
  attr(x, "K1") <- if (length(K1) == 1) K1 else NULL
  attr(x, "delta") <- delta
  attr(x, "dlag") <- dlag
  attr(x, "depvar") <- depvar
  attr(x, "treatvar") <- treatvar
  attr(x, "idvar") <- idvar
  attr(x, "timevar") <- timevar
  class(x) <- c("refbasedmi", class(x))
  x
}

#' Print and summarise RefBasedMI results
#'
#' [RefBasedMI()] returns its stacked long-format imputations as a data frame
#' of class `"refbasedmi"` carrying the run settings as attributes. `print()`
#' shows a one-line description of the run before the data; `summary()`
#' reports the settings together with the row count and the number of missing
#' outcome values in each imputed dataset (which should be zero for every
#' `.imp > 0`).
#'
#' Subsetting keeps the class and the stored settings. If the attributes have
#' been removed, both methods fall back to the plain data-frame behaviour.
#'
#' @param x,object An object returned by [RefBasedMI()].
#' @param ... Passed on to the underlying data-frame method.
#'
#' @return `print()` returns its argument invisibly. `summary()` returns an
#'   object of class `"summary.refbasedmi"`, a list with components `method`,
#'   `reference`, `M`, `K0`, `K1`, `delta`, `dlag`, `depvar`, `n_rows`,
#'   `rows_per_imputation`, and `missing_depvar_per_imputation`.
#'
#' @examples
#' \donttest{
#' asthmaJ2R <- RefBasedMI(data = asthma, depvar = fev, treatvar = treat,
#'                         idvar = id, timevar = time, covar = base,
#'                         method = "J2R", reference = 1, M = 2, seed = 101,
#'                         burnin = 100, verbose = FALSE)
#' summary(asthmaJ2R)
#' head(asthmaJ2R)
#' }
#'
#' @seealso [RefBasedMI()], [as_mids()]
#' @export
print.refbasedmi <- function(x, ...) {
  method <- attr(x, "method")
  if (!is.null(method)) {
    reference <- attr(x, "reference")
    cat("RefBasedMI imputed data: method ", method,
        if (!is.null(reference)) paste0(", reference ", reference),
        ", M = ", attr(x, "M"), "\n", sep = "")
    cat("Long format, ", nrow(x),
        " rows (.imp = 0 holds the original data)\n\n", sep = "")
  }
  NextMethod()
  invisible(x)
}

#' @rdname print.refbasedmi
#' @export
summary.refbasedmi <- function(object, ...) {
  method <- attr(object, "method")
  if (is.null(method)) {
    # attributes were stripped (e.g. by subsetting): plain data-frame summary
    return(NextMethod())
  }
  out <- list(
    method = method,
    reference = attr(object, "reference"),
    M = attr(object, "M"),
    K0 = attr(object, "K0"),
    K1 = attr(object, "K1"),
    delta = attr(object, "delta"),
    dlag = attr(object, "dlag"),
    depvar = attr(object, "depvar"),
    n_rows = nrow(object),
    rows_per_imputation = table(object$.imp),
    missing_depvar_per_imputation =
      tapply(is.na(object[[attr(object, "depvar")]]), object$.imp, sum)
  )
  class(out) <- "summary.refbasedmi"
  out
}

#' @export
print.summary.refbasedmi <- function(x, ...) {
  cat("RefBasedMI imputed data\n")
  cat("  Method:      ", x$method, "\n", sep = "")
  if (!is.null(x$reference)) {
    cat("  Reference:   ", x$reference, "\n", sep = "")
  }
  cat("  Imputations: M = ", x$M, "\n", sep = "")
  if (!is.null(x$K0)) {
    cat("  Causal K0:   ", x$K0, "\n", sep = "")
    cat("  Causal K1:   ", x$K1, "\n", sep = "")
  }
  if (!is.null(x$delta)) {
    cat("  Delta:       ", paste(x$delta, collapse = ", "), "\n", sep = "")
    if (!is.null(x$dlag)) {
      cat("  Delta lag:   ", paste(x$dlag, collapse = ", "), "\n", sep = "")
    }
  }
  cat("  Rows:        ", x$n_rows, " (long format, .imp = 0 is the original data)\n",
      sep = "")
  cat("\nMissing '", x$depvar, "' values by imputation number:\n", sep = "")
  print(x$missing_depvar_per_imputation)
  invisible(x)
}

#' Convert a RefBasedMI result to a mice::mids object
#'
#' A convenience wrapper around [mice::as.mids()] so the imputations can be
#' analysed with [mice::with.mids()] and pooled by Rubin's rules with
#' [mice::pool()]. Requires the suggested package \pkg{mice}.
#'
#' @param x An object returned by [RefBasedMI()] (or any stacked long-format
#'   data frame with an `.imp` column, where `.imp = 0` is the original data).
#'
#' @return A [mice::mids()] object holding the `M` imputed datasets.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("mice", quietly = TRUE)) {
#'   asthmaJ2R <- RefBasedMI(data = asthma, depvar = fev, treatvar = treat,
#'                           idvar = id, timevar = time, covar = base,
#'                           method = "J2R", reference = 1, M = 5, seed = 101,
#'                           burnin = 100, verbose = FALSE)
#'   fit <- with(as_mids(asthmaJ2R),
#'               lm(fev ~ factor(treat), subset = (time == 12)))
#'   summary(mice::pool(fit))
#' }
#' }
#'
#' @seealso [RefBasedMI()], [mice::as.mids()], [mice::pool()]
#' @export
as_mids <- function(x) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("as_mids() requires the 'mice' package; ",
         "install it with install.packages(\"mice\")")
  }
  x <- as.data.frame(x)
  class(x) <- "data.frame"
  mice::as.mids(x)
}
