#' Extract coefficients from a `base_learner`
#'
#' @param object `base_learner` object.
#' @param cause `numeric(1)` or `NULL`. If `NULL`, returns all causes.
#' @param ... Passed to underlying `coef` methods.
#'
#' @return A coefficient object for one cause or a list across causes.
#' @export
coef.base_learner <- function(object, cause= NULL, ...) {

  if (is.null(object$learner_fit)) {
    cat("No fitted model available (learner_fit is NULL).\n")
    return(invisible(object))
  }


  if (is.null(cause)) {
    return(lapply(object$learner_fit, coef, ...))
  } else{
    return(coef(object$learner_fit[[cause]], ...))
  }
}

#' Extract meta-learner coefficients from a fitted ensemble
#'
#' @param object `poisson_superlearner` object.
#' @param cause `numeric(1)` or `NULL`. If `NULL`, returns all causes.
#' @param ... Passed to underlying `coef` methods.
#'
#' @return A coefficient object for one cause or a list across causes.
#' @export
coef.poisson_superlearner <- function(object, cause=NULL,...) {

  if (is.null(object$superlearner)) {
    cat("No fitted model available (learner_fit is NULL).\n")
    return(invisible(object))
  }

  if (is.null(cause)) {

    return(lapply(object$superlearner, function(sl) coef(sl$meta_learner_fit, ...)))

  } else{
    return(coef(object$superlearner[[cause]]$meta_learner_fit, ...))
  }


}

