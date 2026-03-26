#' Extract coefficients from a fitted base learner
#'
#' Convenience method to extract (cause-specific) model coefficients from a fitted
#' `base_learner` returned by [fit_learner()].
#'
#' For competing risks, `fit_learner()` fits one model per cause, stored in
#' `object$learner_fit[[k]]` for `k = 1, 2, ..., K`. This method simply dispatches
#' to the underlying model’s `coef()` method for each fitted object.
#'
#' @param object `base_learner`. A fitted object returned by [fit_learner()].
#' @param cause `numeric(1)` or `NULL`. Which cause to extract coefficients for.
#'   If `NULL`, coefficients are returned for all causes.
#'   Causes are indexed `1, 2, ..., object$data_info$n_crisks` (with `0` reserved for censoring).
#' @param ... Passed to the underlying `coef()` method of the fitted learner object
#'   (learner-dependent; e.g., `s` for `glmnet`).
#'
#' @details
#' **Learner-dependent output.** The returned coefficient object depends on the
#' base learner used (e.g. a numeric vector, a sparse matrix, a list, etc.).
#' This method does not post-process or rename coefficients; it returns the output
#' of `coef(object$learner_fit[[k]], ...)` unchanged.
#'
#' @return
#' If `cause` is a single integer, returns the coefficient object produced by
#' `coef()` for that cause-specific fitted model.
#'
#' If `cause = NULL`, returns a `list` of length `object$data_info$n_crisks`,
#' where element `[[k]]` contains coefficients for cause `k`.
#'
#' If no fitted model is present (`object$learner_fit` is `NULL`), prints a message
#' and returns `invisible(object)`.
#'
#' @examples
#' d <- simulateStenoT1(50, competing_risks = TRUE)
#' lrn <- Learner_glmnet(covariates = c("age", "value_LDL"),
#'                       lambda = 0, cross_validation = FALSE)
#' bl <- fit_learner(d, learner = lrn, id = "id",
#'                   status = "status_cvd", event_time = "time_cvd",
#'                   number_of_nodes = 4)
#'
#' # coefficients for cause 1
#' coef(bl, cause = 1)
#'
#' # coefficients for all causes (list)
#' coef(bl)
#'
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

#' Extract stacking (meta-learner) coefficients from a fitted Poisson Super Learner
#'
#' Extracts the **meta-learner coefficients** (stacking weights) from a fitted
#' `poisson_superlearner` object returned by [Superlearner()].
#'
#' For each cause `k`, the ensemble stores a fitted meta-learner in
#' `object$superlearner[[k]]$meta_learner_fit`. This method dispatches to the
#' underlying `coef()` method for that fitted meta-learner.
#'
#' @param object `poisson_superlearner`. A fitted ensemble returned by [Superlearner()].
#' @param cause `numeric(1)` or `NULL`. Which cause to extract meta-learner
#'   coefficients for. If `NULL`, coefficients are returned for all causes.
#'   Causes are indexed `1, 2, ..., object$data_info$n_crisks`.
#' @param model Scalar model selector. Default is `"sl"` for the stacked super learner.
#'   Other allowed values are:
#'   \describe{
#'     \item{`0` or `"sl"`}{Use the super learner prediction.}
#'     \item{learner label}{Use one stored base learner by its label in
#'       `object$data_info$learners_labels`.}
#'     \item{`"learner_j"`}{Use the `j`-th stored learner.}
#'     \item{integer `j >= 1`}{Use the `j`-th stored learner.}
#'   }
#' @param ... Passed to the underlying `coef()` method of the fitted meta-learner
#'   (learner-dependent; e.g., `s` for `glmnet`).
#'
#' @details
#' **What coefficients represent.** These coefficients correspond to the meta-learner
#' regression of the outcome on the cross-validated base-learner predictions
#' (`Z1`, `Z2`, ...). Under the default meta-learner, they are the stacking
#' weights (on the scale defined by the meta-learner).
#'
#' **Learner-dependent output.** The returned coefficient object depends on the
#' meta-learner implementation (by default a `glmnet` fit, often returning a sparse
#' matrix). This method does not rename `Z*` terms or post-process coefficients; it
#' returns the output of `coef(object$superlearner[[k]]$meta_learner_fit, ...)`
#' unchanged.
#'
#' **Single-learner special case.** If the ensemble was fit with only one base learner,
#' no meta-learner is fit and `meta_learner_fit` is `NULL`. In that case, `coef()`
#' for the `poisson_superlearner` does not have meta-learner coefficients to return.
#'
#' @return
#' If `cause` is a single integer, returns the coefficient object produced by
#' `coef()` for the cause-specific fitted meta-learner.
#'
#' If `cause = NULL`, returns a `list` of length `object$data_info$n_crisks`,
#' where element `[[k]]` contains meta-learner coefficients for cause `k`.
#'
#' If no fitted ensemble is present (`object$superlearner` is `NULL`), prints a message
#' and returns `invisible(object)`.
#'
#' @examples
#' d <- simulateStenoT1(50, competing_risks = TRUE)
#' learners <- list(
#'   glm = Learner_glmnet(covariates = c("age", "value_LDL"), lambda = 0, cross_validation = FALSE),
#'   gam = Learner_gam(covariates = c("age", "value_LDL"))
#' )
#' fit <- Superlearner(d, id="id", status="status_cvd", event_time="time_cvd",
#'                     learners=learners, number_of_nodes=4, nfold=2)
#'
#' # meta-learner coefficients (cause 1)
#' coef(fit, cause = 1)
#'
#' # meta-learner coefficients for all causes (list)
#' coef(fit)
#'
#' @export
coef.poisson_superlearner <- function(object, cause = NULL, model = "sl", ...) {

  if (is.null(object$superlearner)) {
    cat("No fitted model available (superlearner is NULL).\n")
    return(invisible(object))
  }

  if (is.null(cause)) {
    out <- lapply(seq_len(object$data_info$n_crisks), function(k) {
      fit_info <- psl_get_stored_fit(object, cause = k, model = model)
      stats::coef(fit_info$fit, ...)
    })
    names(out) <- paste0("cause_", seq_len(object$data_info$n_crisks))
    return(out)
  }

  fit_info <- psl_get_stored_fit(object, cause = cause, model = model)
  stats::coef(fit_info$fit, ...)
}
