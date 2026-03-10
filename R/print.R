#' Print method for `base_learner`
#'
#' Prints a compact description of the fitted base learner, including the learner
#' type, the time-grid used, and (optionally) the fitted model object for a given
#' cause.
#'
#' @param x `base_learner` object returned by [fit_learner()].
#' @param cause `numeric(1)` or `NULL`. Which cause to print the fitted model for.
#'   If `NULL`, prints one line per cause (classes only) instead of printing the full
#'   fitted objects.
#' @param ... Passed to the underlying fitted object `print()` method when `cause`
#'   is a single integer.
#'
#' @return Invisibly returns `x`.
#' @export
print.base_learner <- function(x, cause=1, ...) {

  if (is.null(x$learner_fit)) {
    cat("No fitted model available (learner_fit is NULL).\n")
    return(invisible(x))
  }

  if (is.null(cause)) {
    return(invisible(lapply(x$learner_fit, print, ...)))
  }

  return(print(x$learner_fit[[cause]], ...))
}


#' Print method for `poisson_superlearner`
#'
#' Prints a compact description of the fitted Poisson Super Learner, including the
#' number of base learners, the meta-learner, the time-grid used, and competing-risk
#' structure. Optionally prints the fitted meta-learner for a given cause.
#'
#' @param x `poisson_superlearner` object returned by [Superlearner()].
#' @param cause `numeric(1)` or `NULL`. Which cause’s meta-learner fit to print.
#'   If `NULL`, prints one line per cause (classes only) instead of printing the full
#'   fitted objects.
#' @param model Scalar model selector. Default is `"sl"` for the stacked super learner.
#'   Other allowed values are:
#'   \describe{
#'     \item{`0` or `"sl"`}{Use the super learner prediction.}
#'     \item{learner label}{Use one stored base learner by its label in
#'       `object$data_info$learners_labels`.}
#'     \item{`"learner_j"`}{Use the `j`-th stored learner.}
#'     \item{integer `j >= 1`}{Use the `j`-th stored learner.}
#'   }
#' @param ... Passed to the underlying fitted meta-learner `print()` method when
#'   `cause` is a single integer.
#'
#' @return Invisibly returns `x`.
#' @export
print.poisson_superlearner <- function(x, cause = 1, model = "sl", ...) {

  if (is.null(x$superlearner)) {
    cat("No fitted model available (superlearner is NULL).\n")
    return(invisible(x))
  }

  if (is.null(cause)) {
    invisible(lapply(seq_len(x$data_info$n_crisks), function(k) {
      fit_info <- psl_get_stored_fit(x, cause = k, model = model)
      print(fit_info$fit, ...)
    }))
    return(invisible(x))
  }

  fit_info <- psl_get_stored_fit(x, cause = cause, model = model)
  print(fit_info$fit, ...)
  invisible(x)
}
