#' Coef method for base_learner objects
#' @export
coef.base_learner <- function(object, cause= NULL, ...) {

  if (is.null(object$learner_fit)) {
    cat("No fitted model available (learner_fit is NULL).\n")
    return(invisible(object))
  }


  if (is.null(cause)) {
    return(lapply(object$learner_fit, coef))
  } else{
    return(coef(object$learner_fit[[cause]], ...))
  }
}

#' Coef method for poisson_superlearner objects
#' @export
coef.poisson_superlearner <- function(object, cause=NULL,...) {

  if (is.null(object$superlearner)) {
    cat("No fitted model available (learner_fit is NULL).\n")
    return(invisible(object))
  }

  if (is.null(cause)) {

    l <- list()
    for(i in length(object$superlearner)){

      l[[i]]<-coef(object$superlearner[[i]]$meta_learner_fit)

    }

    return(l)

  } else{
    return(coef(object$superlearner[[cause]]$meta_learner_fit, ...))
  }


}

