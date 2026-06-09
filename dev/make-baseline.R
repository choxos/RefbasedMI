# Capture reference outputs from the currently loaded RefBasedMI build.
#
# Usage:
#   RBMI_LIB=/path/to/pristine-lib Rscript dev/make-baseline.R [outdir]
#
# Writes one RDS per scenario to outdir (default dev/baseline/) plus a
# sessionInfo.txt describing the environment the baseline was generated in.

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || !nzchar(here)) here <- "dev"
source(file.path(here, "scenarios.R"))

outdir <- if (length(args) >= 1) args[[1]] else file.path(here, "baseline")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Generating baselines into", outdir, "\n")
cat("RefBasedMI version:", as.character(utils::packageVersion("RefBasedMI")), "\n\n")

summary_tbl <- data.frame(scenario = character(), status = character(),
                          stringsAsFactors = FALSE)

for (nm in names(scenarios)) {
  cat(sprintf("  %-26s ... ", nm))
  r <- run_one_scenario(scenarios[[nm]])
  saveRDS(r, file.path(outdir, paste0(nm, ".rds")))
  status <- if (r$ok) "ok" else paste0("ERROR: ", r$error)
  cat(status, "\n")
  summary_tbl <- rbind(summary_tbl,
                       data.frame(scenario = nm, status = status,
                                  stringsAsFactors = FALSE))
}

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
write.csv(summary_tbl, file.path(outdir, "summary.csv"), row.names = FALSE)
cat("\nDone. Wrote", length(scenarios), "baseline files.\n")
