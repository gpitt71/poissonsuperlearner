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
#' If no fitted model is present (`object$learner_fit` is `NULL`), signals a message
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
coef.base_learner <- function(object, cause = NULL, ...) {

  if (is.null(object$learner_fit)) {
    message("No fitted model available (learner_fit is NULL).")
    return(invisible(object))
  }

  if (is.null(cause)) {
    return(lapply(object$learner_fit, coef, ...))
  } else {
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
#' If no fitted ensemble is present (`object$superlearner` is `NULL`), signals a message
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
#' @export
coef.poisson_superlearner <- function(object, cause = NULL, model = "sl", ...) {

  if (is.null(object$superlearner)) {
    message("No fitted model available (superlearner is NULL).")
    return(invisible(object))
  }

  n_crisks <- object$data_info$n_crisks

  learners_by_cause <- object$learners_by_cause
  learners_labels_by_cause <- object$learners_labels_by_cause
  z_covariates_by_cause <- object$z_covariates_by_cause

  ## Backward compatibility with old fitted objects
  if (is.null(learners_by_cause)) {
    learners_by_cause <- replicate(
      n_crisks,
      object$learners,
      simplify = FALSE
    )
  }

  if (is.null(learners_labels_by_cause)) {
    learners_labels_by_cause <- object$data_info$learners_labels_by_cause
  }

  if (is.null(learners_labels_by_cause)) {
    learners_labels_by_cause <- lapply(learners_by_cause, function(x) {
      labs <- names(x)
      if (is.null(labs)) {
        labs <- paste0("learner_", seq_along(x))
      }
      labs
    })
  }

  if (is.null(z_covariates_by_cause)) {
    z_covariates_by_cause <- object$data_info$z_covariates_by_cause
  }

  if (is.null(z_covariates_by_cause)) {
    z_covariates_by_cause <- lapply(
      learners_by_cause,
      function(x) paste0("Z", seq_along(x))
    )
  }

  n_learners_by_cause <- lengths(learners_by_cause)

  if (!is.null(cause)) {
    if (length(cause) != 1L || is.na(cause) || cause != as.integer(cause)) {
      stop("'cause' must be NULL or a single positive integer.", call. = FALSE)
    }

    cause <- as.integer(cause)

    if (cause < 1L || cause > n_crisks) {
      stop(
        sprintf("'cause' must be between 1 and %d.", n_crisks),
        call. = FALSE
      )
    }

    causes_to_extract <- cause
  } else {
    causes_to_extract <- seq_len(n_crisks)
  }


  ## ------------------------------------------------------------
  ## Local model resolver for cause-specific libraries
  ## ------------------------------------------------------------

  resolve_model_for_cause <- function(model, k) {

    if (is.null(model) || length(model) != 1L) {
      stop("'model' must be a scalar.", call. = FALSE)
    }

    if (is.numeric(model)) {
      if (is.na(model) || model != as.integer(model)) {
        stop("Numeric 'model' must be one of 0, 1, 2, ...", call. = FALSE)
      }

      model <- as.integer(model)

      if (model == 0L) {
        return(list(type = "sl", index = 0L, label = "sl"))
      }

      if (model < 1L || model > n_learners_by_cause[k]) {
        stop(
          sprintf(
            "Numeric model %d is unavailable for cause %d. Available learners are 1:%d.",
            model,
            k,
            n_learners_by_cause[k]
          ),
          call. = FALSE
        )
      }

      return(list(
        type = "learner",
        index = model,
        label = learners_labels_by_cause[[k]][model]
      ))
    }

    if (!is.character(model) || is.na(model)) {
      stop("'model' must be a character scalar or a numeric scalar.", call. = FALSE)
    }

    model_chr <- trimws(model)
    model_chr_lc <- tolower(model_chr)

    if (model_chr_lc %in% c("sl", "superlearner", "super_learner")) {
      return(list(type = "sl", index = 0L, label = "sl"))
    }

    if (grepl("^learner_[0-9]+$", model_chr_lc)) {
      j <- as.integer(sub("^learner_", "", model_chr_lc))

      if (j < 1L || j > n_learners_by_cause[k]) {
        stop(
          sprintf(
            "'%s' is unavailable for cause %d. Available learners are learner_1, ..., learner_%d.",
            model_chr,
            k,
            n_learners_by_cause[k]
          ),
          call. = FALSE
        )
      }

      return(list(
        type = "learner",
        index = j,
        label = learners_labels_by_cause[[k]][j]
      ))
    }

    j <- match(model_chr, learners_labels_by_cause[[k]])

    if (is.na(j)) {
      stop(
        sprintf(
          "Learner label '%s' is not available for cause %d.",
          model_chr,
          k
        ),
        call. = FALSE
      )
    }

    list(
      type = "learner",
      index = j,
      label = model_chr
    )
  }


  ## ------------------------------------------------------------
  ## Extract coefficients for one cause
  ## ------------------------------------------------------------

  extract_one_cause <- function(k) {

    sel <- resolve_model_for_cause(model, k)
    sl_k <- object$superlearner[[k]]

    if (sel$type == "sl") {

      fit_k <- sl_k$meta_learner_fit

      if (is.null(fit_k)) {
        return(NULL)
      }

      coefs <- tryCatch(
        psl_extract_meta_coefs(fit_k, ...),
        error = function(e) stats::coef(fit_k, ...)
      )

      zmap <- stats::setNames(
        learners_labels_by_cause[[k]],
        z_covariates_by_cause[[k]]
      )

      if (!is.null(names(coefs))) {
        names(coefs) <- psl_rename_z_in_text(names(coefs), zmap)
      }

      if (is.matrix(coefs) || inherits(coefs, "Matrix")) {
        rn <- rownames(coefs)
        if (!is.null(rn)) {
          rownames(coefs) <- psl_rename_z_in_text(rn, zmap)
        }
      }

      return(coefs)
    }

    j <- sel$index

    fit_j <- if (n_learners_by_cause[k] == 1L) {
      sl_k$learners_fit
    } else {
      sl_k$learners_fit[[j]]
    }

    stats::coef(fit_j, ...)
  }


  ## ------------------------------------------------------------
  ## Return one cause or all causes
  ## ------------------------------------------------------------

  if (length(causes_to_extract) == 1L) {
    return(extract_one_cause(causes_to_extract))
  }

  out <- lapply(causes_to_extract, extract_one_cause)
  names(out) <- paste0("cause_", causes_to_extract)

  out
}
