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
#' @param model Model selector. Default is `"sl"` for the stacked super learner.
#'   Allowed values are:
#'   \describe{
#'     \item{`0`, `"sl"`, `"superlearner"`, or `"super_learner"`}{Summarize the
#'       stacked meta-learner. For causes with no fitted meta-learner, this
#'       falls back to the retained base learner.}
#'     \item{`"discrete_sl"` and aliases}{Summarize the cause-specific base
#'       learners with the smallest cross-validated deviance.}
#'     \item{learner label}{Summarize one stored base learner by its label in
#'       `object$data_info$learners_labels[[k]]`.}
#'     \item{`"learner_j"` or character integer `"j"`}{Summarize the `j`-th
#'       stored learner.}
#'     \item{integer `j >= 1`}{Summarize the `j`-th stored learner.}
#'     \item{vector of labels or positive integer indices}{Use cause-specific
#'       base learners; length must equal `object$data_info$n_crisks`.}
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

  extract_meta_coefs <- function(fit, ...) {

    coefs <- tryCatch(
      stats::coef(fit, s = 0, ...),
      error = function(e1) {
        tryCatch(
          stats::coef(fit, ...),
          error = function(e2) NULL
        )
      }
    )

    if (is.null(coefs)) {
      return(NULL)
    }

    if (is.matrix(coefs) || inherits(coefs, "sparseMatrix")) {
      mat <- as.matrix(coefs)
      out <- as.numeric(mat[, 1L])
      names(out) <- rownames(mat)
      return(out)
    }

    coefs
  }

  rename_z_coefs <- function(coefs, k) {

    if (is.null(coefs) || is.null(names(coefs))) {
      return(coefs)
    }

    labs_k <- learner_labels(k)
    z_names <- paste0("Z", seq_along(labs_k))

    nms <- names(coefs)
    hit <- match(nms, z_names)
    ok <- !is.na(hit)

    nms[ok] <- labs_k[hit[ok]]
    names(coefs) <- nms

    coefs
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

  ## Base learner, discrete SL, or SL fallback without a meta-learner.
  if (model_sel$type == "learner") {

    if (is.null(cause)) {
      out <- lapply(seq_len(n_crisks), function(k) {
        fit_info <- get_stored_fit(k)
        summary(fit_info$fit, ...)
      })
      names(out) <- paste0("cause_", seq_len(n_crisks))
      return(out)
    }

    fit_info <- get_stored_fit(cause)
    return(summary(fit_info$fit, ...))
  }

  causes_to_show <- if (is.null(cause)) {
    seq_len(n_crisks)
  } else {
    cause
  }

  meta_available <- vapply(
    causes_to_show,
    function(k) !is.null(object$superlearner[[k]]$meta_learner_fit),
    logical(1L)
  )

  if (!any(meta_available)) {

    if (is.null(cause)) {
      out <- lapply(seq_len(n_crisks), function(k) {
        fit_info <- get_stored_fit(k)
        summary(fit_info$fit, ...)
      })
      names(out) <- paste0("cause_", seq_len(n_crisks))
      return(out)
    }

    fit_info <- get_stored_fit(cause)
    return(summary(fit_info$fit, ...))
  }

  cat("Call:\n")

  ml_engine <- "<not available>"

  if (!is.null(object$metalearner)) {
    if (!is.null(object$metalearner$engine)) {
      ml_engine <- object$metalearner$engine
    } else {
      ml_engine <- class(object$metalearner)[1L]
    }
  }

  cat(
    sprintf(
      "  Superlearner(..., metalearner = %s)\n",
      ml_engine
    )
  )

  cat("\nFitted object:\n")
  cat("  Class: poisson_superlearner\n")
  cat("  Number of competing risks:", n_crisks, "\n")
  cat("  Number of folds:", object$data_info$nfold, "\n")
  cat("  Maximum follow-up:", object$data_info$maximum_followup, "\n")
  cat("  Number of nodes:", length(object$data_info$nodes), "\n")

  cat("\nRetained learners by cause:\n")
  for (k in seq_len(n_crisks)) {
    labs_k <- learner_labels(k)
    cat(
      sprintf(
        "  cause %d: %s\n",
        k,
        paste(labs_k, collapse = ", ")
      )
    )
  }

  cat("\nCross-validation deviance:\n")
  if (!is.null(object$cross_validation_deviance)) {
    print(object$cross_validation_deviance)
  } else {
    cat("  <not available>\n")
  }


  meta_out <- vector("list", length(causes_to_show))
  names(meta_out) <- paste0("cause_", causes_to_show)

  for (ii in seq_along(causes_to_show)) {

    k <- causes_to_show[ii]
    fit_k <- object$superlearner[[k]]$meta_learner_fit

    if (is.null(fit_k)) {
      fit_info <- get_stored_fit(k)
      cat(
        sprintf(
          "  cause %d: <no meta-learner fitted; using learner '%s'>\n",
          k,
          fit_info$label
        )
      )
      meta_out[[ii]] <- NULL
      next
    }

    coefs <- extract_meta_coefs(fit_k, ...)
    coefs <- rename_z_coefs(coefs, k)

    if (is.null(coefs)) {
      cat("  cause ", k, ": <not available>\n", sep = "")
      meta_out[[ii]] <- NULL
      next
    }

    cat("  cause ", k, ":\n", sep = "")
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
