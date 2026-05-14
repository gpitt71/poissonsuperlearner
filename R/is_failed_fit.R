#' Test whether an object is a failed-fit sentinel
#'
#' Developer note: internal predicate paired with `make_failed_fit()`.
#'
#' @param model Object to test.
#' @returns Logical scalar.
#' @noRd
is_failed_fit <- function(model) {
  inherits(model, "psl_failed_fit")
}
