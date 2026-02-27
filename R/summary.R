#' Summarize a fitted Poisson Super Learner object
#'
#' @param object `poisson_superlearner`.
#' @param ... Unused.
#'
#' @return `data.table` with cross-validated Poisson deviance summaries for
#' learners and meta-learner.
#' @export
summary.poisson_superlearner <- function(object, ...) {

  .get_labels <- function(object) {
    labs <- NULL
    if (!is.null(object$data_info) && !is.null(object$data_info$learners_label)) {
      labs <- object$data_info$learners_label
    }
    if (is.null(labs) || length(labs) == 0L) {
      if (!is.null(object$learners) && length(object$learners) > 0L) {
        labs <- names(object$learners)
      }
    }
    if (is.null(labs) || length(labs) == 0L) {
      labs <- paste0("learner_", seq_along(object$learners))
    }
    labs
  }

  .z_map <- function(labels) {
    stats::setNames(as.character(labels), paste0("Z", seq_along(labels)))
  }

  .rename_z_in_text <- function(x, zmap) {
    if (is.null(x)) return(x)
    x <- as.character(x)
    for (z in names(zmap)) {
      # word boundary avoids changing e.g. Z10 when replacing Z1
      x <- gsub(paste0("\\b", z, "\\b"), zmap[[z]], x, perl = TRUE)
    }
    x
  }

  .extract_meta_coefs <- function(fit) {
    if (is.null(fit)) return(NULL)

    # glmnet / cv.glmnet
    if (inherits(fit, "cv.glmnet")) {
      lam <- fit$lambda.min
      cc <- stats::coef(fit, s = lam)
      cc <- as.matrix(cc)
      out <- cc[, 1]
      return(out)
    }
    if (inherits(fit, "glmnet")) {
      lam <- fit$lambda
      # if multiple lambdas, take the first (or smallest) deterministically
      lam_use <- if (length(lam) >= 1L) lam[[1L]] else NULL
      cc <- if (!is.null(lam_use)) stats::coef(fit, s = lam_use) else stats::coef(fit)
      cc <- as.matrix(cc)
      out <- cc[, 1]
      return(out)
    }

    # fallback
    cc <- tryCatch(stats::coef(fit), error = function(e) NULL)
    if (is.null(cc)) return(NULL)
    cc
  }

  labels <- .get_labels(object)
  zmap <- .z_map(labels)

  cat("Call:\n")
  if (!is.null(object$metalearner)) {
    ml <- object$metalearner

    # Print a compact "call-like" description using fields we know exist
    ml_class <- class(ml)[1]
    labels <- object$data_info$learners_label

    ml_formula <- if (!is.null(labels) && length(labels) > 0L) {
      paste(labels, collapse = " + ")
    } else {
      NULL
    }

    cat("  Meta-learner:\n")
    cat("   ", ml_class, "\n", sep = "")
    if (!is.null(ml_formula)) cat("    ensemble: ", ml_formula, "\n", sep = "")
    if (!is.null(object$data_info$nodes)) cat("    number of time knots: ", length(object$data_info$nodes), "\n", sep = "")
    if (!is.null(object$data_info$nfold)) cat("    number of folds: ", object$data_info$nfold, "\n", sep = "")
    if (!is.null(object$data_info$n_crisks) & object$data_info$n_crisks > 1) cat("    competing risks: ", object$data_info$n_crisks, "\n", sep = "")
    if (!is.null(object$data_info$maximum_followup)) cat("    maximum follow-up: ", object$data_info$maximum_followup, "\n", sep = "")
  } else {
    cat("  Meta-learner: <none>\n")
  }

  cat("\nCross-validation deviance (Average across V-Folds):\n")
  if (!is.null(object$cross_validation_deviance)) {
    print(object$cross_validation_deviance)
  } else if (!is.null(object$meta_learner_cross_validation)) {
    # backward-compatible fallback if your object still uses the older name
    print(object$meta_learner_cross_validation)
  } else {
    cat("  <not available>\n")
  }

  cat("\nMeta-learner coefficients:\n")
  if (is.null(object$superlearner) || length(object$superlearner) == 0L) {
    cat("  <not available>\n")
    return(invisible(object))
  }

  # one set of coefficients per competing risk (k)
  for (k in seq_along(object$superlearner)) {
    fit_k <- object$superlearner[[k]]$meta_learner_fit
    coefs <- .extract_meta_coefs(fit_k)

    if (is.null(coefs)) {
      cat("  k = ", k, ": <not available>\n", sep = "")
      next
    }

    nms <- names(coefs)
    nms2 <- .rename_z_in_text(nms, zmap)
    names(coefs) <- nms2

    cat("  k = ", k, ":\n", sep = "")
    print(coefs)
    cat("\n")
  }

  invisible(object)
}



#' Summary method for base_learner objects
#' Summarize a fitted base learner object
#'
#' @param object `base_learner`.
#' @param cause `numeric(1)`. Cause index.
#' @param ... Unused.
#'
#' @return A list containing model-specific summary information for the selected
#' cause.
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
