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

  pred <- predict(
    object,
    newdata = newdata,
    times = times,
    cause = cause,
    model = model
  )

  if (is.null(pred)) {
    return(NULL)
  }

  time_col <- object$data_info$event_time
  id_col <- object$data_info$id

  for (ix in seq_along(times)) {
    tmp <- pred[get(time_col) == times[ix]]

    if (nrow(tmp) == 0L) {
      next
    }

    data.table::setorderv(tmp, id_col)
    rr_output[, ix] <- tmp[["absolute_risk"]]
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
predictRisk.base_learner <- function(object,
                                     newdata,
                                     times,
                                     cause = 1,
                                     ...) {

  rr_output <- matrix(
    NA_real_,
    nrow = nrow(newdata),
    ncol = length(times)
  )

  pred <- predict(
    object,
    newdata = newdata,
    times = times,
    cause = cause
  )

  if (is.null(pred)) {
    return(NULL)
  }

  time_col <- object$data_info$event_time
  id_col <- object$data_info$id

  for (ix in seq_along(times)) {
    tmp <- pred[get(time_col) == times[ix]]

    if (nrow(tmp) == 0L) {
      next
    }

    data.table::setorderv(tmp, id_col)
    rr_output[, ix] <- tmp[["absolute_risk"]]
  }

  return(rr_output)
}
