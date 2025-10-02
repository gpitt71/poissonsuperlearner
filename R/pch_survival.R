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
#' @examples
#' # Example: 2 subjects, 2 causes
#' id <- c(1,1,2,2)
#' dt <- c(1,2,1,2)
#' haz <- matrix(c(0.1,0.05, 0.2,0.1), ncol = 2, byrow = TRUE)
#' pch_survival(id, dt, haz)
#'
#' @export
#' @useDynLib tmlensemble, .registration = TRUE
#' @importFrom Rcpp evalCpp
pch_survival <- function(id, dt, haz, na_is_zero = FALSE) {
  .Call(`_tmlensemble_pch_survival`, id, dt, haz, na_is_zero)
}
