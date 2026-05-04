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

  n_crisks <- object$data_info$n_crisks

  ## ------------------------------------------------------------
  ## Recover cause-specific bookkeeping, with backward compatibility
  ## ------------------------------------------------------------

  learners_by_cause <- object$learners_by_cause
  learners_labels_by_cause <- object$learners_labels_by_cause
  z_covariates_by_cause <- object$z_covariates_by_cause

  if (is.null(learners_by_cause)) {
    if (is.null(object$learners)) {
      stop("No learner library found in the fitted object.", call. = FALSE)
    }

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

  if (length(learners_by_cause) != n_crisks) {
    stop("'object$learners_by_cause' must have length equal to n_crisks.", call. = FALSE)
  }

  if (length(learners_labels_by_cause) != n_crisks) {
    stop("'object$learners_labels_by_cause' must have length equal to n_crisks.", call. = FALSE)
  }

  if (length(z_covariates_by_cause) != n_crisks) {
    stop("'object$z_covariates_by_cause' must have length equal to n_crisks.", call. = FALSE)
  }


  ## ------------------------------------------------------------
  ## Validate cause
  ## ------------------------------------------------------------

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

    causes_to_show <- cause
  } else {
    causes_to_show <- seq_len(n_crisks)
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
  ## If model selects a base learner, dispatch to that fitted learner
  ## ------------------------------------------------------------

  first_sel <- resolve_model_for_cause(model, causes_to_show[1L])

  if (first_sel$type == "learner") {

    out <- lapply(causes_to_show, function(k) {
      sel_k <- resolve_model_for_cause(model, k)
      sl_k <- object$superlearner[[k]]

      fit_k <- if (n_learners_by_cause[k] == 1L) {
        sl_k$learners_fit
      } else {
        sl_k$learners_fit[[sel_k$index]]
      }

      summary(fit_k, ...)
    })

    names(out) <- paste0("cause_", causes_to_show)

    if (length(out) == 1L) {
      return(out[[1L]])
    }

    return(out)
  }


  ## ------------------------------------------------------------
  ## Default: summarize the stacked superlearner
  ## ------------------------------------------------------------

  cat("Call:\n")
  cat("  Superlearner(...)\n")

  cat("\nFitted object:\n")
  cat("  Class: poisson_superlearner\n")
  cat("  Number of competing risks:", n_crisks, "\n")
  cat("  Number of folds:", object$data_info$nfold, "\n")
  cat("  Maximum follow-up:", object$data_info$maximum_followup, "\n")
  cat("  Number of nodes:", length(object$data_info$nodes), "\n")

  cat("\nLearners by cause:\n")
  for (k in seq_len(n_crisks)) {
    cat(
      sprintf(
        "  cause %d: %d learner(s): %s\n",
        k,
        n_learners_by_cause[k],
        paste(learners_labels_by_cause[[k]], collapse = ", ")
      )
    )
  }

  cat("\nMeta-learner by cause:\n")
  for (k in seq_len(n_crisks)) {
    fit_k <- object$superlearner[[k]]$meta_learner_fit

    if (is.null(fit_k)) {
      cat(
        sprintf(
          "  cause %d: <none; direct learner: %s>\n",
          k,
          learners_labels_by_cause[[k]][1L]
        )
      )
    } else {
      cat(
        sprintf(
          "  cause %d: %s\n",
          k,
          class(fit_k)[1L]
        )
      )
    }
  }

  cat("\nCross-validation deviance (average across V-folds):\n")
  if (!is.null(object$cross_validation_deviance)) {
    print(object$cross_validation_deviance)
  } else if (!is.null(object$meta_learner_cross_validation)) {
    print(object$meta_learner_cross_validation)
  } else {
    cat("  <not available>\n")
  }


  ## ------------------------------------------------------------
  ## Meta-learner coefficients
  ## ------------------------------------------------------------

  cat("\nMeta-learner coefficients:\n")

  meta_out <- vector("list", length(causes_to_show))
  names(meta_out) <- paste0("cause_", causes_to_show)

  for (ii in seq_along(causes_to_show)) {

    k <- causes_to_show[ii]
    fit_k <- object$superlearner[[k]]$meta_learner_fit

    if (is.null(fit_k)) {
      cat(
        sprintf(
          "  cause %d: <not available; direct learner: %s>\n",
          k,
          learners_labels_by_cause[[k]][1L]
        )
      )
      meta_out[[ii]] <- NULL
      next
    }

    coefs <- tryCatch(
      psl_extract_meta_coefs(fit_k, ...),
      error = function(e) NULL
    )

    if (is.null(coefs)) {
      cat(sprintf("  cause %d: <not available>\n", k))
      meta_out[[ii]] <- NULL
      next
    }

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

    cat(sprintf("  cause %d:\n", k))
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
