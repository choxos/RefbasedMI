# Internal replacement for Hmisc::summarize(X, by = llist(...), FUN = sum).
#
# Returns a data frame with one row per unique combination of the grouping
# variables (the grouping columns first, then the summarised value), ordered by
# the grouping variables in the order supplied, with the first varying slowest.
# This row order is load-bearing: callers take cumsum() over the result to build
# the cumulative case counts that index the imputation loops, so the ordering
# must match Hmisc::summarize() exactly.

group_summarize <- function(x, by, FUN = sum) {
  if (!is.list(by)) {
    stop("'by' must be a list of grouping variables")
  }
  if (is.null(names(by)) || any(!nzchar(names(by)))) {
    names(by) <- paste0("by", seq_along(by))
  }
  agg <- stats::aggregate(list(value = x), by = by, FUN = FUN)
  agg <- agg[do.call(order, agg[names(by)]), , drop = FALSE]
  rownames(agg) <- NULL
  agg
}
