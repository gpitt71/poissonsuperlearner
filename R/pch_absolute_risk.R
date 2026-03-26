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
#' id <- c(1L, 1L, 2L, 2L)
#' dt <- c(1, 1, 1, 1)
#' haz <- rbind(
#'   c(0.10, 0.05),
#'   c(0.20, 0.10),
#'   c(0.05, 0.02),
#'   c(0.10, 0.03)
#' )
#' pch_absolute_risk(id = id, dt = dt, haz = haz, cause_idx = 1)
#'
#'
#' @export
#' @useDynLib poissonsuperlearner, .registration = TRUE
#' @importFrom Rcpp evalCpp
pch_absolute_risk <- function(id, dt, haz, cause_idx, one_based = TRUE, na_is_zero = FALSE) {
  .Call(`_poissonsuperlearner_pch_absolute_risk`, id, dt, haz, as.integer(cause_idx),
        one_based, na_is_zero)
}


#' Absolute risk (Euler approximation) for a cause under piecewise-constant hazards
#'
#' Computes the cumulative incidence function using the first-order Euler
#' (discrete) approximation:
#'   \deqn{F_j(t) \approx \sum S(t_{k-1}) \lambda_{j,k} \Delta t_k}{
#'     F_j(t) = approx sum S(t_{k-1}) lambda_{j,k} Delta t_k
#'   }
#' Grouped by `id`, this returns the cumulative incidence at the end of each interval.
#'
#' @param id Integer vector. Sorted by `id` then time.
#' @param dt Numeric vector of interval lengths.
#' @param haz Numeric matrix (n x C) of cause-specific hazards per interval.
#' @param cause_idx Integer. Index of the cause of interest (1-based by default).
#' @param one_based Logical. If `TRUE`, `cause_idx` is 1-based. If `FALSE`, 0-based.
#' @param na_is_zero Logical. If `TRUE`, treat NA/Inf hazards as zero.
#'
#' @return Numeric vector of cumulative incidence values (Euler approximation)
#'         at the end of each interval.
#'
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
#' pch_absolute_risk_euler(id = id, dt = dt, haz = haz, cause_idx = 1)
#'
#' @export
#' @useDynLib poissonsuperlearner, .registration = TRUE
#' @importFrom Rcpp evalCpp
pch_absolute_risk_euler <- function(id, dt, haz, cause_idx, one_based = TRUE, na_is_zero = FALSE) {
  .Call(`_poissonsuperlearner_pch_absolute_risk_euler`, id, dt, haz, as.integer(cause_idx),
        one_based, na_is_zero)
}










