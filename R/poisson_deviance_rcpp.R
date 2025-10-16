#' @title Poisson deviance for piece-wise constant exponential models
#' @description
#' Computes the Poisson deviance for survival or competing-risk models
#' with one or more hazards, using the compiled C++ routine
#' \code{poisson_deviance_cpp()}.
#'
#' @param dt_learners A \code{data.table} containing at least the columns:
#'   \itemize{
#'     \item \code{folder} – integer or factor grouping for cross-validation folds.
#'     \item \code{learner} – character identifier of the model.
#*     \item \code{tij} – numeric time interval width for each observation.
#'     \item \code{pwch_1}, …, \code{pwch_H} – hazard rates.
#'     \item \code{delta_1}, …, \code{delta_H} – event indicators.
#'   }
#' @param pwch_cols Character vector of hazard column names.
#' @param delta_cols Character vector of corresponding event indicator columns.
#'
#' @return A \code{data.table} with columns \code{learner} and mean \code{deviance}.
#'
#' @export
#' @useDynLib tmlensemble, .registration = TRUE
#' @importFrom Rcpp evalCpp
poisson_deviance_rcpp <- function(dt_learners, pwch_cols, delta_cols) {
  stopifnot(data.table::is.data.table(dt_learners))
  H <- length(pwch_cols)
  if (length(delta_cols) != H) stop("delta_cols must match pwch_cols length")

  lambda <- as.matrix(dt_learners[, ..pwch_cols])
  delta  <- as.matrix(dt_learners[, ..delta_cols])
  storage.mode(delta) <- "integer"

  res <- poisson_deviance_cpp(
    folder  = dt_learners[["folder"]],
    learner = dt_learners[["learner"]],
    lambda  = lambda,
    tij     = dt_learners[["tij"]],
    delta   = delta
  )

  data.table::as.data.table(res)[order(learner)]
}
