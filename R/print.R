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

  ## ------------------------------------------------------------
  ## Recover cause-specific bookkeeping, with backward compatibility
  ## ------------------------------------------------------------

  learners_by_cause <- x$learners_by_cause
  learners_labels_by_cause <- x$learners_labels_by_cause
  z_covariates_by_cause <- x$z_covariates_by_cause

  if (is.null(learners_by_cause)) {
    if (is.null(x$learners)) {
      stop("No learner library found in the fitted object.", call. = FALSE)
    }

    learners_by_cause <- replicate(
      n_crisks,
      x$learners,
      simplify = FALSE
    )
  }

  if (is.null(learners_labels_by_cause)) {
    learners_labels_by_cause <- x$data_info$learners_labels_by_cause
  }

  if (is.null(learners_labels_by_cause)) {
    learners_labels_by_cause <- lapply(learners_by_cause, function(z) {
      labs <- names(z)
      if (is.null(labs)) {
        labs <- paste0("learner_", seq_along(z))
      }
      labs
    })
  }

  if (is.null(z_covariates_by_cause)) {
    z_covariates_by_cause <- x$data_info$z_covariates_by_cause
  }

  if (is.null(z_covariates_by_cause)) {
    z_covariates_by_cause <- lapply(
      learners_by_cause,
      function(z) paste0("Z", seq_along(z))
    )
  }

  n_learners_by_cause <- lengths(learners_by_cause)

  if (length(learners_by_cause) != n_crisks) {
    stop("'x$learners_by_cause' must have length equal to n_crisks.", call. = FALSE)
  }

  if (length(learners_labels_by_cause) != n_crisks) {
    stop("'x$learners_labels_by_cause' must have length equal to n_crisks.", call. = FALSE)
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

    causes_to_print <- cause
  } else {
    causes_to_print <- seq_len(n_crisks)
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
  ## Retrieve the stored fit for one cause
  ## ------------------------------------------------------------

  get_fit_for_cause <- function(k) {

    sel <- resolve_model_for_cause(model, k)
    sl_k <- x$superlearner[[k]]

    if (sel$type == "sl") {

      if (!is.null(sl_k$meta_learner_fit)) {
        return(list(
          fit = sl_k$meta_learner_fit,
          type = "sl",
          label = "sl"
        ))
      }

      ## Direct learner fallback when no meta-learner exists for this cause
      return(list(
        fit = sl_k$learners_fit,
        type = "learner",
        label = learners_labels_by_cause[[k]][1L]
      ))
    }

    j <- sel$index

    fit_j <- if (n_learners_by_cause[k] == 1L) {
      sl_k$learners_fit
    } else {
      sl_k$learners_fit[[j]]
    }

    list(
      fit = fit_j,
      type = "learner",
      label = sel$label
    )
  }


  ## ------------------------------------------------------------
  ## Compact print for all causes
  ## ------------------------------------------------------------

  if (is.null(cause)) {

    cat("Poisson Super Learner\n")
    cat("  Number of competing risks:", n_crisks, "\n")
    cat("  Number of folds:", x$data_info$nfold, "\n")
    cat("  Maximum follow-up:", x$data_info$maximum_followup, "\n")
    cat("  Number of nodes:", length(x$data_info$nodes), "\n")

    cat("\nStored fits:\n")

    for (k in causes_to_print) {
      fit_info <- get_fit_for_cause(k)

      fit_class <- if (is.null(fit_info$fit)) {
        "<NULL>"
      } else {
        paste(class(fit_info$fit), collapse = ", ")
      }

      cat(
        sprintf(
          "  cause %d: model = %s, label = %s, class = %s\n",
          k,
          fit_info$type,
          fit_info$label,
          fit_class
        )
      )
    }

    return(invisible(x))
  }


  ## ------------------------------------------------------------
  ## Full print for one selected cause
  ## ------------------------------------------------------------

  fit_info <- get_fit_for_cause(cause)

  cat("Poisson Super Learner\n")
  cat("  Cause:", cause, "\n")
  cat("  Model:", fit_info$type, "\n")
  cat("  Label:", fit_info$label, "\n\n")

  print(fit_info$fit, ...)

  invisible(x)
}
