# Shared regression scenario definitions for RefBasedMI.
#
# Both dev/make-baseline.R (captures outputs from the pristine v0.2.0 build) and
# dev/compare-baseline.R (re-runs against the working tree) source this file so
# the two halves can never drift apart. Each scenario is a zero-argument closure
# that calls RefBasedMI() with literal bare names, because the function captures
# depvar/treatvar/idvar/timevar/covar/method/methodvar via substitute().
#
# The library to load RefBasedMI from is taken from the RBMI_LIB environment
# variable (falls back to the default .libPaths()).

rbmi_lib <- Sys.getenv("RBMI_LIB", unset = NA)
if (!is.na(rbmi_lib) && nzchar(rbmi_lib)) {
  library(RefBasedMI, lib.loc = rbmi_lib)
} else {
  library(RefBasedMI)
}

# Bundled datasets (deterministic, identical across both runs).
data(asthma, package = "RefBasedMI")
data(antidepressant, package = "RefBasedMI")
data(acupuncture, package = "RefBasedMI")

# S14 fixture: punch interim (non-monotone) holes into asthma. For a fixed set of
# ids that are observed at the final visit (time == 12), blank the time == 8 value
# so it becomes an interim missing value feeding the fillinterims() path.
asthma_interim <- asthma
.interim_ids <- sort(unique(asthma_interim$id))[1:15]
asthma_interim$fev[asthma_interim$id %in% .interim_ids & asthma_interim$time == 8] <- NA

scenarios <- list(
  S01_asthma_MAR = function() RefBasedMI(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = "MAR", M = 2, seed = 101, burnin = 100),

  S02_asthma_J2R = function() RefBasedMI(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = "J2R", reference = 1, M = 2, seed = 101, burnin = 100),

  S03_asthma_CR = function() RefBasedMI(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = "CR", reference = 1, M = 2, seed = 101, burnin = 100),

  S04_asthma_CIR = function() RefBasedMI(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = "CIR", reference = 1, M = 2, seed = 101, burnin = 100),

  S05_asthma_LMCF = function() RefBasedMI(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = "LMCF", M = 2, seed = 101, burnin = 100),

  S06_asthma_Causal = function() RefBasedMI(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = "Causal", reference = 1, K0 = 1, K1 = 0.5,
    M = 2, seed = 101, burnin = 100),

  S07_asthma_J2R_ref2 = function() RefBasedMI(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = "J2R", reference = 2, M = 2, seed = 101, burnin = 100),

  S08_asthma_J2R_delta = function() RefBasedMI(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = "J2R", reference = 1, M = 2, seed = 101, burnin = 100,
    delta = c(0.5, 0.5, 1, 1), dlag = c(1, 1, 1, 1)),

  S09_asthma_J2R_nocovar = function() RefBasedMI(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    method = "J2R", reference = 1, M = 2, seed = 101, burnin = 100),

  S10_antidep_methodvar = function() RefBasedMI(
    data = antidepressant, depvar = HAMD17.TOTAL, treatvar = TREATMENT.NAME,
    idvar = PATIENT.NUMBER, timevar = VISIT.NUMBER, covar = basval,
    methodvar = methodcol, referencevar = referencecol,
    M = 2, seed = 101, burnin = 100),

  S11_acupuncture_LMCF = function() RefBasedMI(
    data = acupuncture, depvar = head, treatvar = treat, idvar = id, timevar = time,
    covar = head_base, method = "LMCF", M = 2, seed = 101, burnin = 100),

  S12_acupuncture_J2R = function() RefBasedMI(
    data = acupuncture, depvar = head, treatvar = treat, idvar = id, timevar = time,
    covar = head_base, method = "J2R", reference = 1, M = 2, seed = 101, burnin = 100),

  S13_asthma_J2R_mle = function() RefBasedMI(
    data = asthma, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = "J2R", reference = 1, M = 2, seed = 101, burnin = 100,
    mle = TRUE),

  S14_asthma_interim = function() RefBasedMI(
    data = asthma_interim, depvar = fev, treatvar = treat, idvar = id, timevar = time,
    covar = base, method = "J2R", reference = 1, M = 2, seed = 101, burnin = 100)
)

# Run one scenario quietly; return list(ok, value|error).
run_one_scenario <- function(fn) {
  res <- tryCatch(
    suppressWarnings(suppressMessages(
      utils::capture.output(out <- fn())
    )),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    list(ok = FALSE, error = conditionMessage(res))
  } else {
    list(ok = TRUE, value = out)
  }
}
