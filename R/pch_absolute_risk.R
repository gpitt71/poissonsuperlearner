#' Absolute risk (cumulative incidence) for a cause under piecewise-constant hazards
#'
#' Computes, per row, the cumulative incidence function at the end of each interval,
#' grouped by `id`. The number of causes is inferred from the number of columns in `haz`.
#'
#' @param id Integer vector. Sorted by `id` then time.
#' @param dt Numeric vector of interval lengths.
#' @param haz Numeric matrix (n x C) of cause-specific hazards per interval.
#'           Columns correspond to causes 1..C.
#' @param cause_idx Integer. Index of the cause of interest (1-based by default).
#' @param one_based Logical. If `TRUE`, `cause_idx` is 1-based. If `FALSE`, 0-based.
#' @param na_is_zero Logical. If `TRUE`, treat NA/Inf hazards as zero.
#'
#' @return Numeric vector of cumulative incidence values at the end of each interval.
#'
#' @examples
#' # Hazard columns pwch_1, pwch_2, ...
#' cols <- grep("^pwch_[0-9]+$", names(dt_long), value = TRUE)
#' haz  <- as.matrix(dt_long[, ..cols])
#' # cause 2:
#' cif2 <- pch_absolute_risk(dt_long$id, dt_long$deltatime, haz, cause_idx = 2)
#'
#' @export
#' @useDynLib tmlensemble, .registration = TRUE
#' @importFrom Rcpp evalCpp
pch_absolute_risk <- function(id, dt, haz, cause_idx, one_based = TRUE, na_is_zero = FALSE) {
  .Call(`_tmlensemble_pch_absolute_risk`, id, dt, haz, as.integer(cause_idx),
        one_based, na_is_zero)
}
