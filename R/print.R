#' Print method for base_learner
#' @export
print.base_learner <- function(x, cause=1, ...) {

  if (is.null(x$learner_fit)) {
    cat("No fitted model available (learner_fit is NULL).\n")
    return(invisible(x))
  }

  return(print(x$learner_fit[[cause]], ...))
}

#' Print method for poisson_superlearner objects
#' @export
print.poisson_superlearner <- function(object, cause=1, ...) {

  if (is.null(object$superlearner)) {
    cat("No fitted model available (learner_fit is NULL).\n")
    return(invisible(object))
  }


  return(print(object$superlearner[[cause]]$meta_learner_fit, ...))
}
