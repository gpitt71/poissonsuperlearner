#' Print method for `base_learner`
#'
#' Prints a compact description of the fitted base learner, including the learner
#' type, the time-grid used, and (optionally) the fitted model object for a given
#' cause.
#'
#' @param x `base_learner` object returned by [fit_learner()].
#' @param cause `numeric(1)` or `NULL`. Which cause to print the fitted model for.
#'   If `NULL`, prints one line per cause (classes only) instead of printing the full
#'   fitted objects.
#' @param ... Passed to the underlying fitted object `print()` method when `cause`
#'   is a single integer.
#'
#' @return Invisibly returns `x`.
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
#' print(bl, cause = NULL)
#'
#'
#' @export
print.base_learner <- function(x, cause=1, ...) {

  if (is.null(x$learner_fit)) {
    cat("No fitted model available (learner_fit is NULL).\n")
    return(invisible(x))
  }

  if (is.null(cause)) {
    return(invisible(lapply(x$learner_fit, print, ...)))
  }

  return(print(x$learner_fit[[cause]], ...))
}


#' Print method for `poisson_superlearner`
#'
#' Prints a compact description of the fitted Poisson Super Learner, including the
#' number of base learners, the meta-learner, the time-grid used, and competing-risk
#' structure. Optionally prints the fitted meta-learner for a given cause.
#'
#' @param x `poisson_superlearner` object returned by [Superlearner()].
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
#' @param ... Passed to the underlying fitted meta-learner `print()` method when
#'   `cause` is a single integer.
#'
#' @return Invisibly returns `x`.
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
#' print(fit, cause = NULL)
#'
#' @export
print.poisson_superlearner <- function(x, cause = 1, model = "sl", ...) {

  if (is.null(x$superlearner) || length(x$superlearner) == 0L) {
    cat("No fitted model available (superlearner is NULL).\n")
    return(invisible(x))
  }

  n_crisks <- x$data_info$n_crisks
  learners_by_cause <- x$learners

  if (
    is.null(learners_by_cause) ||
    !is.list(learners_by_cause) ||
    length(learners_by_cause) != n_crisks
  ) {
    stop("Could not find cause-specific learners in `x$learners`.", call. = FALSE)
  }

  learner_labels <- function(k) {

    labs <- NULL

    if (
      !is.null(x$data_info$learners_labels) &&
      is.list(x$data_info$learners_labels) &&
      length(x$data_info$learners_labels) >= k
    ) {
      labs <- x$data_info$learners_labels[[k]]
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

    cv_dt <- x$cross_validation_deviance

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
          "`model = 'discrete_sl'` requires `x$cross_validation_deviance`, unless each cause has only one retained learner.",
          call. = FALSE
        )
      }

      cv_df <- as.data.frame(cv_dt)

      if (!("deviance" %in% names(cv_df))) {
        stop(
          "`x$cross_validation_deviance` must contain a `deviance` column for `model = 'discrete_sl'`.",
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
          "`x$cross_validation_deviance` must contain `cause_index` or cause labels matching `x$learners` for `model = 'discrete_sl'`.",
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

    sl_k <- x$superlearner[[k]]
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

  if (is.null(cause)) {
    invisible(lapply(seq_len(n_crisks), function(k) {
      fit_info <- get_stored_fit(k)
      cat("\nCause ", k, ", model = ", fit_info$label, "\n", sep = "")
      print(fit_info$fit, ...)
    }))
    return(invisible(x))
  }

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

  fit_info <- get_stored_fit(as.integer(cause))
  print(fit_info$fit, ...)

  invisible(x)
}
