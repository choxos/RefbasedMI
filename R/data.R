#' Asthma trial data
#'
#' The asthma clinical trial data used in the Stata \emph{mimix} help file. The
#' outcome \code{fev} (forced expiratory volume) is recorded in long format at
#' weeks 2, 4, 8 and 12, with some post-baseline values missing.
#'
#' @format A data frame with 732 rows and 5 variables:
#'   \describe{
#'   \item{id}{patient identifier}
#'   \item{time}{measurement time point, in weeks (2, 4, 8, 12)}
#'   \item{treat}{randomised treatment arm (1 or 2)}
#'   \item{base}{baseline value of the outcome, used as a covariate}
#'   \item{fev}{forced expiratory volume, the outcome (missing after discontinuation)}
#'   }
#' @source Distributed with the Stata \emph{mimix} package by Suzie Cro.
"asthma"


#' Antidepressant trial data
#'
#' The antidepressant clinical trial data described by White, Royes and Best
#' (2020). The outcome \code{HAMD17.TOTAL} (Hamilton depression rating) is
#' recorded in long format at visits 4, 5, 6 and 7. The data also carry
#' per-patient method and reference columns used to illustrate individual-specific
#' imputation.
#'
#' @format A data frame with 688 rows and 14 variables:
#'   \describe{
#'   \item{PATIENT.NUMBER}{patient identifier}
#'   \item{HAMA.TOTAL}{Hamilton anxiety rating scale total}
#'   \item{PGI.IMPROVEMENT}{patient global impression of improvement}
#'   \item{VISIT.DATE}{day of the visit relative to randomisation}
#'   \item{VISIT.NUMBER}{visit number (4, 5, 6, 7)}
#'   \item{TREATMENT.NAME}{randomised treatment arm}
#'   \item{PATIENT.SEX}{patient sex}
#'   \item{POOLED.INVESTIGATOR}{pooled investigator centre}
#'   \item{basval}{baseline value of the outcome, used as a covariate}
#'   \item{HAMD17.TOTAL}{Hamilton 17-item depression rating, the outcome}
#'   \item{change}{change from baseline in the outcome}
#'   \item{miss_flag}{indicator that the record is observed}
#'   \item{methodcol}{per-patient imputation method (for individual-specific imputation)}
#'   \item{referencecol}{per-patient reference arm (for individual-specific imputation)}
#'   }
#' @source White I R, Royes J, Best N (2020) \doi{10.1080/10543406.2019.1684308}.
"antidepressant"


#' Acupuncture trial data
#'
#' Data from a randomised, double-blind, parallel-group trial comparing active
#' treatment with placebo. The outcome \code{head} (headache severity score) is
#' recorded in long format at times 3 and 12.
#'
#' @format A data frame with 802 rows and 11 variables:
#'   \describe{
#'   \item{id}{patient identifier}
#'   \item{time}{measurement time point (3 or 12)}
#'   \item{age}{patient age, in years}
#'   \item{sex}{patient sex}
#'   \item{migraine}{indicator of migraine (versus tension-type headache)}
#'   \item{chronicity}{duration of headache disorder, in years}
#'   \item{practice_id}{recruiting practice identifier}
#'   \item{treat}{randomised treatment arm (1 or 2)}
#'   \item{head_base}{baseline headache score, used as a covariate}
#'   \item{head}{headache severity score, the outcome (missing after discontinuation)}
#'   \item{withdrawal_reason}{reason for withdrawal, where recorded}
#'   }
#' @source Vickers A J et al. acupuncture for chronic headache trial.
"acupuncture"
