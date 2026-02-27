#' Print method for `base_learner`
#'
#' @param x `base_learner` object.
#' @param cause `numeric(1)`. Cause index.
#' @param ... Unused.
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
#' @param object `poisson_superlearner` object.
#' @param cause `numeric(1)`. Cause index.
#' @param ... Unused.
#' @return Invisibly returns `object`.
#' @export
print.poisson_superlearner <- function(object, cause=1, ...) {

  if (is.null(object$superlearner)) {
    cat("No fitted model available (learner_fit is NULL).\n")
    return(invisible(object))
  }


  if (is.null(cause)) {
    return(invisible(lapply(object$superlearner, function(sl) print(sl$meta_learner_fit, ...))))
  }

  return(print(object$superlearner[[cause]]$meta_learner_fit, ...))
}

