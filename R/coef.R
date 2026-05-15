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
#' @param model Model selector. Default is `"sl"` for the stacked super learner.
#'   Allowed values are:
#'   \describe{
#'     \item{`0`, `"sl"`, `"superlearner"`, or `"super_learner"`}{Extract
#'       coefficients from the stacked meta-learner. For causes with no fitted
#'       meta-learner, this falls back to the retained base learner.}
#'     \item{`"discrete_sl"` and aliases}{Extract coefficients from the
#'       cause-specific base learners with the smallest cross-validated
#'       deviance.}
#'     \item{learner label}{Extract coefficients from one stored base learner by
#'       its label in `object$data_info$learners_labels[[k]]`.}
#'     \item{`"learner_j"` or character integer `"j"`}{Extract coefficients from
#'       the `j`-th stored learner.}
#'     \item{integer `j >= 1`}{Extract coefficients from the `j`-th stored
#'       learner.}
#'     \item{vector of labels or positive integer indices}{Use cause-specific
#'       base learners; length must equal `object$data_info$n_crisks`.}
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
#' `coef()` for the selected cause-specific fitted model: the meta-learner when
#' `model = "sl"` and a meta-learner is available, or the selected base learner
#' when `model` selects a base learner or no meta-learner is available.
#'
#' If `cause = NULL`, returns a `list` of length `object$data_info$n_crisks`,
#' where element `[[k]]` contains coefficients for the selected model for cause
#' `k`.
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
coef.poisson_superlearner <- function(object, cause = NULL, model = "sl", ...) {

  if (is.null(object$superlearner) || length(object$superlearner) == 0L) {
    message("No fitted model available (superlearner is NULL).")
    return(invisible(object))
  }

  n_crisks <- object$data_info$n_crisks
  learners_by_cause <- object$learners

  if (
    is.null(learners_by_cause) ||
    !is.list(learners_by_cause) ||
    length(learners_by_cause) != n_crisks
  ) {
    stop("Could not find cause-specific learners in `object$learners`.", call. = FALSE)
  }

  learner_labels <- function(k) {

    labs <- NULL

    if (
      !is.null(object$data_info$learners_labels) &&
      is.list(object$data_info$learners_labels) &&
      length(object$data_info$learners_labels) >= k
    ) {
      labs <- object$data_info$learners_labels[[k]]
    }

    if (
      is.null(labs) ||
      length(labs) != length(learners_by_cause[[k]]) ||
      anyNA(labs) ||
      any(!nzchar(labs))
    ) {
      labs <- names(learners_by_cause[[k]])
    }

    if (
      is.null(labs) ||
      length(labs) != length(learners_by_cause[[k]]) ||
      anyNA(labs) ||
      any(!nzchar(labs))
    ) {
      labs <- paste0("learner_", seq_along(learners_by_cause[[k]]))
    }

    labs
  }

  best_discrete_sl_index_by_cause <- function() {

    cv_dt <- object$cross_validation_deviance

    index_by_cause <- integer(n_crisks)
    label_by_cause <- character(n_crisks)

    cause_names <- names(learners_by_cause)

    for (k in seq_len(n_crisks)) {

      labs_k <- learner_labels(k)
      n_k <- length(labs_k)

      if (n_k == 1L) {
        index_by_cause[k] <- 1L
        label_by_cause[k] <- labs_k[1L]
        next
      }

      if (
        is.null(cv_dt) ||
        !is.data.frame(cv_dt) ||
        nrow(cv_dt) == 0L
      ) {
        stop(
          "`model = 'discrete_sl'` requires `object$cross_validation_deviance`, unless each cause has only one retained learner.",
          call. = FALSE
        )
      }

      cv_df <- as.data.frame(cv_dt)

      if (!("deviance" %in% names(cv_df))) {
        stop(
          "`object$cross_validation_deviance` must contain a `deviance` column for `model = 'discrete_sl'`.",
          call. = FALSE
        )
      }

      if ("cause_index" %in% names(cv_df)) {
        cv_k <- cv_df[as.integer(cv_df[["cause_index"]]) == k, , drop = FALSE]
      } else if (
        "cause" %in% names(cv_df) &&
        !is.null(cause_names) &&
        length(cause_names) >= k &&
        !is.na(cause_names[k]) &&
        nzchar(cause_names[k])
      ) {
        cv_k <- cv_df[as.character(cv_df[["cause"]]) == cause_names[k], , drop = FALSE]
      } else {
        stop(
          "`object$cross_validation_deviance` must contain `cause_index` or cause labels matching `object$learners` for `model = 'discrete_sl'`.",
          call. = FALSE
        )
      }

      if (nrow(cv_k) == 0L) {
        stop(
          sprintf(
            "No cross-validation deviance is available for cause %d, so `model = 'discrete_sl'` cannot choose among %d learners.",
            k,
            n_k
          ),
          call. = FALSE
        )
      }

      dev_k <- as.numeric(cv_k[["deviance"]])

      if (all(is.na(dev_k) | !is.finite(dev_k))) {
        stop(
          sprintf(
            "All cross-validation deviances are missing or non-finite for cause %d.",
            k
          ),
          call. = FALSE
        )
      }

      dev_k[is.na(dev_k) | !is.finite(dev_k)] <- Inf
      best_row <- which.min(dev_k)

      ix <- NA_integer_

      if ("learner" %in% names(cv_k)) {
        best_label <- as.character(cv_k[["learner"]][best_row])
        ix <- match(best_label, labs_k)
      }

      if (is.na(ix) && "learner_index" %in% names(cv_k)) {
        ix <- as.integer(cv_k[["learner_index"]][best_row])
      }

      if (is.na(ix) || ix < 1L || ix > n_k) {
        stop(
          sprintf(
            "The best cross-validation learner for cause %d could not be matched to the retained learner library. Available learners: %s.",
            k,
            paste(sprintf("%d='%s'", seq_along(labs_k), labs_k), collapse = ", ")
          ),
          call. = FALSE
        )
      }

      index_by_cause[k] <- ix
      label_by_cause[k] <- labs_k[ix]
    }

    list(
      type = "learner",
      index_by_cause = index_by_cause,
      label_by_cause = label_by_cause
    )
  }

  resolve_model <- function(model) {

    if (is.null(model) || !length(model)) {
      stop(
        "`model` must be 'sl', 'discrete_sl', a learner label, a learner index, or a vector of cause-specific learner labels/indices.",
        call. = FALSE
      )
    }

    if (!(length(model) %in% c(1L, n_crisks))) {
      stop(
        sprintf(
          "`model` must have length 1 or length %d, one entry per cause.",
          n_crisks
        ),
        call. = FALSE
      )
    }

    if (is.numeric(model) || is.integer(model)) {

      if (anyNA(model) || any(model != as.integer(model))) {
        stop("Numeric `model` values must be integer learner indices.", call. = FALSE)
      }

      model_int <- as.integer(model)

      if (length(model_int) == 1L && model_int == 0L) {
        return(list(
          type = "sl",
          index_by_cause = rep(0L, n_crisks),
          label_by_cause = rep("sl", n_crisks)
        ))
      }

      if (any(model_int < 1L)) {
        stop(
          "Numeric `model` values must be positive learner indices. Use scalar 0 or 'sl' for the Super Learner.",
          call. = FALSE
        )
      }

      index_by_cause <- if (length(model_int) == 1L) {
        rep(model_int, n_crisks)
      } else {
        model_int
      }

      for (k in seq_len(n_crisks)) {
        labs_k <- learner_labels(k)

        if (index_by_cause[k] > length(labs_k)) {
          stop(
            sprintf(
              "Learner index %d is not available for cause %d. Available learners: %s.",
              index_by_cause[k],
              k,
              paste(sprintf("%d='%s'", seq_along(labs_k), labs_k), collapse = ", ")
            ),
            call. = FALSE
          )
        }
      }

      return(list(
        type = "learner",
        index_by_cause = index_by_cause,
        label_by_cause = vapply(
          seq_len(n_crisks),
          function(k) learner_labels(k)[index_by_cause[k]],
          character(1L)
        )
      ))
    }

    if (!is.character(model)) {
      stop(
        "`model` must be 'sl', 'discrete_sl', a learner label, a learner index, or a vector of cause-specific learner labels/indices.",
        call. = FALSE
      )
    }

    if (anyNA(model) || any(!nzchar(trimws(model)))) {
      stop("Character `model` values cannot be missing or empty.", call. = FALSE)
    }

    model_chr <- trimws(model)
    model_lc <- tolower(model_chr)

    sl_values <- c("sl", "superlearner", "super_learner")

    discrete_sl_values <- c(
      "discrete_sl",
      "discrete_superlearner",
      "discrete_super_learner",
      "discretesl"
    )

    if (length(model_chr) == 1L && model_lc %in% sl_values) {
      return(list(
        type = "sl",
        index_by_cause = rep(0L, n_crisks),
        label_by_cause = rep("sl", n_crisks)
      ))
    }

    if (length(model_chr) == 1L && model_lc %in% discrete_sl_values) {
      return(best_discrete_sl_index_by_cause())
    }

    if (any(model_lc %in% c(sl_values, discrete_sl_values))) {
      stop(
        "`model = 'sl'` and `model = 'discrete_sl'` must be supplied as scalar model selectors. Do not mix them with cause-specific base-learner selectors.",
        call. = FALSE
      )
    }

    model_chr_by_cause <- if (length(model_chr) == 1L) {
      rep(model_chr, n_crisks)
    } else {
      model_chr
    }

    index_by_cause <- integer(n_crisks)
    label_by_cause <- character(n_crisks)

    for (k in seq_len(n_crisks)) {

      labs_k <- learner_labels(k)
      selector_k <- model_chr_by_cause[k]
      selector_k_lc <- tolower(selector_k)

      if (grepl("^[0-9]+$", selector_k)) {
        ix <- as.integer(selector_k)
      } else if (grepl("^learner_[0-9]+$", selector_k_lc)) {
        ix <- as.integer(sub("^learner_", "", selector_k_lc))
      } else {
        ix <- match(selector_k, labs_k)
      }

      if (is.na(ix) || ix < 1L || ix > length(labs_k)) {
        stop(
          sprintf(
            "Learner '%s' is not available for cause %d. Available learners: %s.",
            selector_k,
            k,
            paste(sprintf("%d='%s'", seq_along(labs_k), labs_k), collapse = ", ")
          ),
          call. = FALSE
        )
      }

      index_by_cause[k] <- ix
      label_by_cause[k] <- labs_k[ix]
    }

    list(
      type = "learner",
      index_by_cause = index_by_cause,
      label_by_cause = label_by_cause
    )
  }

  model_sel <- resolve_model(model)

  get_stored_fit <- function(k) {

    sl_k <- object$superlearner[[k]]
    n_k <- length(learners_by_cause[[k]])

    if (
      model_sel$type == "sl" &&
      !is.null(sl_k$meta_learner_fit)
    ) {
      return(list(
        fit = sl_k$meta_learner_fit,
        type = "sl",
        index = 0L,
        label = "sl"
      ))
    }

    learner_ix <- if (model_sel$type == "learner") {
      model_sel$index_by_cause[k]
    } else {
      1L
    }

    fits_k <- sl_k$learners_fit

    fit_k <- if (n_k == 1L) {
      fits_k
    } else {
      fits_k[[learner_ix]]
    }

    list(
      fit = fit_k,
      type = "learner",
      index = learner_ix,
      label = learner_labels(k)[learner_ix]
    )
  }

  if (!is.null(cause)) {
    if (
      length(cause) != 1L ||
      is.na(cause) ||
      cause != as.integer(cause) ||
      cause < 1L ||
      cause > n_crisks
    ) {
      stop(
        sprintf("`cause` must be NULL or a single integer between 1 and %d.", n_crisks),
        call. = FALSE
      )
    }

    cause <- as.integer(cause)
  }

  if (is.null(cause)) {
    out <- lapply(seq_len(n_crisks), function(k) {
      fit_info <- get_stored_fit(k)
      stats::coef(fit_info$fit, ...)
    })
    names(out) <- paste0("cause_", seq_len(n_crisks))
    return(out)
  }

  fit_info <- get_stored_fit(cause)
  stats::coef(fit_info$fit, ...)
}
