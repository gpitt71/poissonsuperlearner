#' Create a lightweight failed-fit sentinel
#'
#' Developer note: internal helper used to avoid retaining heavy failed model
#' objects in learner libraries.
#'
#' The returned object records a failure reason and has classes used by
#' `is_failed_fit()`. It has no side effects.
#'
#' @param reason Character description of the failure.
#' @returns A `psl_failed_fit` sentinel object.
#' @noRd
make_failed_fit <- function(reason = NA_character_) {
  structure(
    list(reason = reason),
    class = c("psl_failed_fit", "psl_fit_stub")
  )
}
