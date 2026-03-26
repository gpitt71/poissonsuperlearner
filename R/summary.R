#' Summarize a fitted Poisson Super Learner object
#'
#' Prints:
#' 1) a compact description of the fitted ensemble,
#' 2) cross-validated deviances for base learners (when available),
#' 3) cause-specific meta-learner coefficients (stacking weights).
#'
#' @param object `poisson_superlearner` returned by [Superlearner()].
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
#' @param ... Passed to the underlying `coef()` method for the fitted meta-learner
#'   (learner-dependent; e.g. `s` for `glmnet`).
#'
#' @return Invisibly returns a `list` with elements:
#' \describe{
#'   \item{cross_validation_deviance}{`data.table` (or `NULL`).}
#'   \item{meta_coefficients}{List of length `n_crisks` with cause-specific coefficient objects (or `NULL`).}
#' }
#'
#' @examples
#' d <- simulateStenoT1(30, competing_risks = TRUE)
#'
#' learners <- list(
#'   lasso = Learner_glmnet(
#'     covariates = "sex",
#'     alpha = 1,
#'     lambda = 0.01,
#'     cross_validation = FALSE
#'   ),
#'   ridge = Learner_glmnet(
#'     covariates = c("sex", "value_LDL"),
#'     alpha = 0,
#'     lambda = 0.01,
#'     cross_validation = FALSE
#'   )
#' )
#'
#' fit <- Superlearner(
#'   data = d,
#'   id = "id",
#'   status = "status_cvd",
#'   event_time = "time_cvd",
#'   learners = learners,
#'   number_of_nodes = 3,
#'   nfold = 2
#' )
#'
#' s <- summary(fit, cause = 1)
#' names(s)
#'
#'
#' @export
summary.poisson_superlearner <- function(object,
                                         cause = NULL,
                                         model = "sl",
                                         ...) {

  if (is.null(object$superlearner) || length(object$superlearner) == 0L) {
    cat("No fitted model available (superlearner is NULL).\n")
    return(invisible(object))
  }

  model_sel <- resolve_prediction_model(object, model)
  labels <- psl_get_labels(object)

  ## If the user asks for a stored base learner, or if there is no meta-learner
  ## (single-learner special case), dispatch directly to the underlying fit.
  if (model_sel$type == "learner" || is.null(object$metalearner)) {

    if (is.null(cause)) {
      out <- lapply(seq_len(object$data_info$n_crisks), function(k) {
        fit_info <- psl_get_stored_fit(object, cause = k, model = model)
        summary(fit_info$fit, ...)
      })
      names(out) <- paste0("cause_", seq_len(object$data_info$n_crisks))
      return(out)
    }

    fit_info <- psl_get_stored_fit(object, cause = cause, model = model)
    return(summary(fit_info$fit, ...))
  }

  ## Default: summarize the stacked superlearner
  zmap <- stats::setNames(labels, paste0("Z", seq_along(labels)))

  cat("Call:\n")
  if (!is.null(object$metalearner)) {
    ml <- object$metalearner
    ml_class <- class(ml)[1]
    cat(
      sprintf(
        "  Superlearner(..., learners = c(%s), metalearner = %s)\n",
        paste(labels, collapse = ", "),
        ml_class
      )
    )
  } else {
    cat("  Superlearner(...)\n")
  }

  cat("\nFitted object:\n")
  cat("  Class: poisson_superlearner\n")
  cat("  Number of competing risks:", object$data_info$n_crisks, "\n")
  cat("  Number of learners:", length(object$learners), "\n")
  cat("  Learners:", paste(labels, collapse = ", "), "\n")
  cat("  Number of folds:", object$data_info$nfold, "\n")
  cat("  Maximum follow-up:", object$data_info$maximum_followup, "\n")
  cat("  Number of nodes:", length(object$data_info$nodes), "\n")

  cat("\nCross-validation deviance (Average across V-Folds):\n")
  if (!is.null(object$cross_validation_deviance)) {
    print(object$cross_validation_deviance)
  } else if (!is.null(object$meta_learner_cross_validation)) {
    print(object$meta_learner_cross_validation)
  } else {
    cat("  <not available>\n")
  }

  cat("\nMeta-learner coefficients:\n")

  causes_to_show <- if (is.null(cause)) {
    seq_along(object$superlearner)
  } else {
    cause
  }

  meta_out <- vector("list", length(causes_to_show))
  names(meta_out) <- paste0("cause_", causes_to_show)

  for (ii in seq_along(causes_to_show)) {
    k <- causes_to_show[ii]
    fit_k <- object$superlearner[[k]]$meta_learner_fit

    coefs <- tryCatch(
      psl_extract_meta_coefs(fit_k, ...),
      error = function(e) NULL
    )

    if (is.null(coefs)) {
      cat("  k = ", k, ": <not available>\n", sep = "")
      meta_out[[ii]] <- NULL
      next
    }

    nms <- names(coefs)
    names(coefs) <- psl_rename_z_in_text(nms, zmap)

    cat("  k = ", k, ":\n", sep = "")
    print(coefs)
    cat("\n")

    meta_out[[ii]] <- coefs
  }

  invisible(list(
    cross_validation_deviance = object$cross_validation_deviance,
    meta_coefficients = meta_out
  ))
}


#' Summarize a fitted base learner object
#'
#' Dispatches to the underlying fitted model’s `summary()` method for the selected
#' cause, or returns a list of summaries for all causes.
#'
#' @param object `base_learner` returned by [fit_learner()].
#' @param cause `numeric(1)` or `NULL`. Which cause to summarize. If `NULL`,
#'   returns one summary per cause.
#' @param ... Passed to the underlying `summary()` method (learner-dependent).
#'
#' @return If `cause` is a single integer, returns the underlying model summary for
#' that cause. If `cause = NULL`, returns a list of summaries (one per cause).
#'
#' @examples
#' d <- simulateStenoT1(30, competing_risks = TRUE)
#' lrn <- Learner_glmnet(
#'   covariates = c("sex", "value_LDL"),
#'   lambda = 0.01,
#'   cross_validation = FALSE
#' )
#' bl <- fit_learner(
#'   d,
#'   learner = lrn,
#'   id = "id",
#'   status = "status_cvd",
#'   event_time = "time_cvd",
#'   number_of_nodes = 3
#' )
#' out <- summary(bl, cause = 1)
#'
#'
#' @export
summary.base_learner <- function(object, cause=1, ...) {

  if (is.null(object$learner_fit)) {
    cat("No fitted model available (learner_fit is NULL).\n")
    return(invisible(object))
  }

  if (is.null(cause)) {
    return(lapply(object$learner_fit, summary, ...))
  } else{

    return(summary(object$learner_fit[[cause]], ...))

  }

}
