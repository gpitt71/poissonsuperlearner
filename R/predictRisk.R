#' Absolute-risk matrix predictions for a fitted Poisson Super Learner
#'
#' S3 method compatible with `riskRegression::predictRisk` returning one column
#' per requested time.
#'
#' @param object `poisson_superlearner`. Fitted object.
#' @param newdata `data.frame`. New covariate data.
#' @param times `numeric`. Prediction times.
#' @param cause `numeric(1)`. Cause index.
#' @param model Scalar model selector. Default is `"sl"`. Allowed values are the
#' same as in [predict.poisson_superlearner()].
#' @param ... Unused.
#'
#' @return `numeric` matrix with `nrow(newdata)` rows and `length(times)` columns.
#' @export
#' @export predictRisk.poisson_superlearner
predictRisk.poisson_superlearner <- function(object,
                                             newdata,
                                             times,
                                             cause = 1,
                                             model = "sl",
                                             ...) {

  rr_output <- matrix(
    NA_real_,
    nrow = nrow(newdata),
    ncol = length(times)
  )

  for (ix in seq_along(times)) {
    rr_output[, ix] <- predict(
      object,
      newdata = newdata,
      times = times[ix],
      cause = cause,
      model = model
    )$absolute_risk
  }

  return(rr_output)
}

#' Absolute-risk matrix predictions for a fitted base learner
#'
#' @param object `base_learner`. Fitted object from [fit_learner()].
#' @param newdata `data.frame`. New covariate data.
#' @param times `numeric`. Prediction times.
#' @param cause `numeric(1)`. Cause index.
#' @param ... Unused.
#'
#' @return `numeric` matrix with `nrow(newdata)` rows and `length(times)` columns.
#' @export
predictRisk.base_learner <- function(object,newdata,times,cause=1, ...){

  rr_output <- matrix(NA_real_,
                      nrow=nrow(newdata),
                      ncol=length(times))


  for(ix in seq_along(times)){


    rr_output[,ix]<- predict(object,
                             newdata = newdata,
                             times = times[ix],
                             cause = cause)$absolute_risk

  }

  return(rr_output)


}
