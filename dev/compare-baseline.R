# Re-run every scenario against the working-tree build and diff against the
# captured baselines. Exits non-zero if any scenario diverges, so it can gate
# each refactoring phase.
#
# Usage:
#   RBMI_LIB=/path/to/dev-lib Rscript dev/compare-baseline.R [baselinedir]

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || !nzchar(here)) here <- "dev"
source(file.path(here, "scenarios.R"))

basedir <- if (length(args) >= 1) args[[1]] else file.path(here, "baseline")
tol <- as.numeric(Sys.getenv("RBMI_TOL", unset = "1e-12"))

cat("Comparing working tree against baselines in", basedir,
    "(tolerance", tol, ")\n")
cat("RefBasedMI version:", as.character(utils::packageVersion("RefBasedMI")), "\n\n")

n_pass <- 0L; n_fail <- 0L; n_skip <- 0L
failures <- character()

for (nm in names(scenarios)) {
  bf <- file.path(basedir, paste0(nm, ".rds"))
  cat(sprintf("  %-26s ... ", nm))
  if (!file.exists(bf)) { cat("SKIP (no baseline)\n"); n_skip <- n_skip + 1L; next }
  base <- readRDS(bf)
  now <- run_one_scenario(scenarios[[nm]])

  if (!base$ok && !now$ok) {
    # Both error: require the same message (a fixed bug may legitimately change this).
    if (identical(base$error, now$error)) { cat("ok (both error, same msg)\n"); n_pass <- n_pass + 1L }
    else { cat("CHANGED ERROR\n    was:", base$error, "\n    now:", now$error, "\n")
           n_fail <- n_fail + 1L; failures <- c(failures, nm) }
    next
  }
  if (base$ok != now$ok) {
    cat("STATUS FLIP (baseline ok=", base$ok, ", now ok=", now$ok, ")\n", sep = "")
    if (!now$ok) cat("    now error:", now$error, "\n")
    n_fail <- n_fail + 1L; failures <- c(failures, nm); next
  }

  # 0.3.1 classes the result and carries run settings as attributes; the gate
  # compares the data itself, so normalise both sides to a bare data frame
  strip <- function(v) {
    if (is.data.frame(v)) {
      v <- as.data.frame(v)
      attributes(v) <- attributes(v)[c("names", "row.names", "class")]
      class(v) <- "data.frame"
    }
    v
  }
  cmp <- all.equal(strip(base$value), strip(now$value), tolerance = tol)
  if (isTRUE(cmp)) { cat("identical\n"); n_pass <- n_pass + 1L }
  else { cat("DIFF\n"); for (l in cmp) cat("    ", l, "\n")
         n_fail <- n_fail + 1L; failures <- c(failures, nm) }
}

cat(sprintf("\n%d passed, %d failed, %d skipped\n", n_pass, n_fail, n_skip))
if (n_fail > 0) {
  cat("FAILED scenarios:", paste(failures, collapse = ", "), "\n")
  quit(status = 1L, save = "no")
}
cat("All scenarios match baseline.\n")
