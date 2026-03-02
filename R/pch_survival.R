#' Piecewise-constant hazards survival function
#'
#' Computes survival at the end of each interval for competing risks with
#' piecewise constant hazards.
#'
#' @param id Integer vector of subject IDs, sorted by id then time.
#' @param dt Numeric vector of interval lengths.
#' @param haz Numeric matrix (n x C) of cause-specific hazards.
#' @param na_is_zero Logical. If TRUE, treat NA hazards as zero.
#'
#' @return Numeric vector of survival probabilities at the end of each interval.
#'
#' @export
#' @useDynLib poissonsuperlearner, .registration = TRUE
#' @importFrom Rcpp evalCpp
pch_survival <- function(id, dt, haz, na_is_zero = FALSE) {
  .Call(`_poissonsuperlearner_pch_survival`, id, dt, haz, na_is_zero)
}
