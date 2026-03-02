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
#' @param ... Passed to the underlying fitted meta-learner `print()` method when
#'   `cause` is a single integer.
#'
#' @return Invisibly returns `x`.
#' @export
print.poisson_superlearner <- function(x, cause=1, ...) {

  if (is.null(x$superlearner)) {
    cat("No fitted model available (learner_fit is NULL).\n")
    return(invisible(x))
  }


  if (is.null(cause)) {
    return(invisible(lapply(x$superlearner, function(sl) print(sl$meta_learner_fit, ...))))
  }

  return(print(x$superlearner[[cause]]$meta_learner_fit, ...))
}

