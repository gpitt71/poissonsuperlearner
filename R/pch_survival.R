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
#' id <- c(1L, 1L, 2L, 2L)
#' dt <- c(1, 1, 1, 1)
#' haz <- rbind(
#'   c(0.10, 0.05),
#'   c(0.20, 0.10),
#'   c(0.05, 0.02),
#'   c(0.10, 0.03)
#' )
#' pch_survival(id = id, dt = dt, haz = haz)
#'
#' @export
#' @useDynLib poissonsuperlearner, .registration = TRUE
#' @importFrom Rcpp evalCpp
pch_survival <- function(id, dt, haz, na_is_zero = FALSE) {
  .Call(`_poissonsuperlearner_pch_survival`, id, dt, haz, na_is_zero)
}
