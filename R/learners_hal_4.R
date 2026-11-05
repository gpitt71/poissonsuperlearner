#' Token-only HAL Poisson learner with automatic lambda paths
#'
#' `Learner_hal4` is the cross-validated counterpart of [Learner_hal_3()]. It
#' uses the same tokenized HAL basis and Poisson IRLS coordinate-descent solver,
#' but generates a glmnet-style lambda path when `lambda` is not supplied.
#' Cross-validation is always performed.
#'
#' @inheritParams Learner_hal_3
#' @param nlambda Number of values in an automatically generated lambda path.
#' @param lambda_min_ratio Ratio of the smallest to the largest automatically
#'   generated lambda. If `NULL`, uses `1e-4` when the number of observations is
#'   at least the number of penalized basis functions and `1e-2` otherwise.
#' @param lambda Optional explicit decreasing lambda sequence. When supplied,
#'   automatic sequence generation is skipped and this path is used by all CV
#'   folds.
#' @param ... Arguments passed to the tokenized Poisson solver. An explicit
#'   decreasing `lambda` sequence is supported and is shared across CV folds.
#' @return An initialized `Learner_hal4` Reference Class object.
#' @export Learner_hal4
#' @exportClass Learner_hal4
Learner_hal4 <- methods::setRefClass(
  "Learner_hal4", contains = "Learner_hal_3", methods = list(
    initialize = function(covariates = NA_character_, intercept = TRUE,
                          max_degree = 2L, num_knots = c(10L, 5L),
                          maxit_prefit = NA_real_, standardize = TRUE,
                          nlambda = 100L, lambda_min_ratio = NULL,
                          lambda = NULL, ...) {
      callSuper(covariates = covariates, intercept = intercept,
                cross_validation = FALSE, max_degree = max_degree,
                num_knots = num_knots, maxit_prefit = maxit_prefit,
                standardize = standardize, lambda = lambda, ...)
      .self$cross_validation <- TRUE
      .self$fit_arguments$nlambda <- nlambda
      if (!is.null(lambda_min_ratio))
        .self$fit_arguments$lambda_min_ratio <- lambda_min_ratio
    },

    private_fit = function(data, cause, grid_nodes, ...) {
      .self$private_fit_all_causes(
        data, causes = cause, grid_nodes = grid_nodes, ...)[[1L]]
    },

    private_fit_all_causes = function(data, causes, grid_nodes, ...) {
      tryCatch({
        prep <- hal2_prepare_training_data(
          data, .self$basic_covariates, causes, grid_nodes)
        tokens <- hal3_tokens(prep$data, c(.self$basic_covariates, "node"),
                              .self$max_degree, .self$num_knots)
        args <- .self$fit_arguments
        supported <- c(
          "lambda", "nlambda", "lambda_min_ratio", "nfolds", "foldid",
          "tol", "maxit_irls", "maxit_cd", "active_set", "strong_rules",
          "rescale_penalty_factor", "standardize"
        )
        unknown <- setdiff(names(args), c(supported, "intercept"))
        if (length(unknown))
          stop("Unsupported options: ", paste(unknown, collapse = ", "))
        if (is.null(args$nfolds)) args$nfolds <- 5L

        ans <- setNames(vector("list", length(causes)), as.character(causes))
        for (i in seq_along(causes)) {
          fit <- do.call(cv_hal4_fit_tokens, c(list(
            tokens = tokens, data = prep$data, y = prep$y[, i],
            offset = log(prep$data$tij),
            penalty_factor = tokens$default_penalty_factor,
            intercept = .self$intercept
          ), args[intersect(names(args), supported)]))
          .self$lambda_opt <- fit$lambda.min
          ans[[i]] <- fit
        }
        ans
      }, error = function(e) {
        setNames(rep(list(make_failed_fit(conditionMessage(e))), length(causes)),
                 as.character(causes))
      })
    }
  )
)

# Calculate the exact penalty factors used by hal3_fit_tokens(). Keeping this
# local to hal4 avoids changing the established Learner_hal_3 implementation.
hal4_effective_penalty <- function(tokens, penalty_factor, intercept,
                                   rescale_penalty_factor, standardize) {
  p <- tokens$ncol
  penalty_factor <- as.numeric(penalty_factor)
  if (length(penalty_factor) != p || any(!is.finite(penalty_factor)) ||
      any(penalty_factor < 0))
    stop("Invalid `penalty_factor`.")
  if (!any(penalty_factor > 0)) stop("All penalties are zero.")
  penalty_scale <- if (isTRUE(rescale_penalty_factor))
    p / sum(penalty_factor) else 1
  column_scale <- rep(1, p)
  if (isTRUE(standardize)) {
    prevalence <- hal_token_sums_cached_cpp(
      tokens$index, tokens$numeric_data, tokens$cat_data,
      rep(1, nrow(tokens$numeric_data)), rep(0, nrow(tokens$numeric_data)))
    prevalence <- as.numeric(prevalence$curvature)
    column_scale <- if (isTRUE(intercept))
      sqrt(pmax(0, prevalence - prevalence^2)) else sqrt(prevalence)
  }
  penalty_factor * penalty_scale * column_scale
}

hal4_lambda_sequence <- function(tokens, y, offset,
                                 penalty_factor = NULL, intercept = TRUE,
                                 tol = 1e-8, maxit_irls = 100L,
                                 maxit_cd = 10000L, active_set = TRUE,
                                 strong_rules = TRUE,
                                 rescale_penalty_factor = TRUE,
                                 standardize = FALSE, nlambda = 100L,
                                 lambda_min_ratio = NULL) {
  n <- nrow(tokens$numeric_data)
  p <- tokens$ncol
  if (is.null(penalty_factor)) penalty_factor <- rep(1, p)
  effective_pf <- hal4_effective_penalty(
    tokens, penalty_factor, intercept, rescale_penalty_factor, standardize)
  penalized <- is.finite(effective_pf) & effective_pf > 0
  p_penalized <- sum(penalized)
  if (!p_penalized) stop("No penalized candidate basis functions.")

  if (length(nlambda) != 1L || is.na(nlambda) || !is.finite(nlambda) ||
      nlambda < 1 || nlambda != as.integer(nlambda))
    stop("`nlambda` must be one positive integer.", call. = FALSE)
  nlambda <- as.integer(nlambda)
  if (is.null(lambda_min_ratio))
    lambda_min_ratio <- if (n >= p_penalized) 1e-4 else 1e-2
  if (length(lambda_min_ratio) != 1L || is.na(lambda_min_ratio) ||
      !is.finite(lambda_min_ratio) || lambda_min_ratio <= 0 ||
      lambda_min_ratio > 1)
    stop("`lambda_min_ratio` must be in (0, 1].", call. = FALSE)

  # The pilot invokes the unchanged solver so its intercept and zero-penalty
  # coefficients are the actual null model for this objective and offset.
  pilot <- hal3_fit_tokens(
    tokens, y, offset, lambda = 1e30, penalty_factor = penalty_factor,
    intercept = intercept, tol = tol, maxit_irls = maxit_irls,
    maxit_cd = maxit_cd, active_set = active_set,
    strong_rules = strong_rules,
    rescale_penalty_factor = rescale_penalty_factor,
    standardize = standardize)
  score <- hal_token_score_cpp(
    tokens$index, tokens$numeric_data, tokens$cat_data, as.numeric(y),
    as.numeric(offset), pilot$beta[, 1L], pilot$a0[1L])
  lambda_max <- max(abs(score[penalized]) / effective_pf[penalized])
  if (!is.finite(lambda_max))
    stop("Could not calculate a finite `lambda_max`.", call. = FALSE)
  if (lambda_max <= sqrt(.Machine$double.eps))
    return(rep(0, nlambda))
  if (nlambda == 1L) return(lambda_max)
  exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
          length.out = nlambda))
}

hal4_fit_tokens <- function(tokens, y, offset, lambda = NULL,
                            penalty_factor = NULL, intercept = TRUE,
                            tol = 1e-8, maxit_irls = 100L,
                            maxit_cd = 10000L, active_set = TRUE,
                            strong_rules = TRUE,
                            rescale_penalty_factor = TRUE,
                            standardize = FALSE, nlambda = 100L,
                            lambda_min_ratio = NULL) {
  if (is.null(lambda) || !length(lambda)) {
    lambda <- hal4_lambda_sequence(
      tokens, y, offset, penalty_factor, intercept, tol, maxit_irls,
      maxit_cd, active_set, strong_rules, rescale_penalty_factor,
      standardize, nlambda, lambda_min_ratio)
  }
  hal3_fit_tokens(
    tokens, y, offset, lambda = lambda, penalty_factor = penalty_factor,
    intercept = intercept, tol = tol, maxit_irls = maxit_irls,
    maxit_cd = maxit_cd, active_set = active_set,
    strong_rules = strong_rules,
    rescale_penalty_factor = rescale_penalty_factor,
    standardize = standardize)
}

hal4_subset_tokens <- function(tokens, rows) {
  out <- tokens
  out$numeric_data <- tokens$numeric_data[rows, , drop = FALSE]
  out$cat_data <- tokens$cat_data[rows, , drop = FALSE]
  out$index <- hal_token_index_cpp(
    out$numeric_data, out$cat_data, out$token_offsets, out$token_var,
    out$token_kind, out$token_cut, out$token_level)
  out$native_index_bytes <- hal_token_index_bytes_cpp(out$index)
  out
}

cv_hal4_fit_tokens <- function(tokens, data, y, offset, lambda = NULL,
                               penalty_factor = NULL, intercept = TRUE,
                               nfolds = 5L, foldid = NULL, ...) {
  n <- nrow(tokens$numeric_data)
  if (nrow(data) != n || length(y) != n || length(offset) != n)
    stop("Mismatched cross-validation data dimensions.")
  if (length(nfolds) != 1L || is.na(nfolds) || nfolds < 2L ||
      nfolds > n || nfolds != as.integer(nfolds))
    stop("`nfolds` must be an integer between 2 and the sample size.")
  nfolds <- as.integer(nfolds)
  if (is.null(foldid)) foldid <- sample(rep(seq_len(nfolds), length.out = n))
  foldid <- as.integer(foldid)
  if (length(foldid) != n || !setequal(unique(foldid), seq_len(nfolds)))
    stop("`foldid` must assign every row to one of all requested folds.")

  full <- hal4_fit_tokens(
    tokens, y, offset, lambda = lambda, penalty_factor = penalty_factor,
    intercept = intercept, ...)
  cvraw <- matrix(NA_real_, nfolds, length(full$lambda))
  for (fold in seq_len(nfolds)) {
    training <- foldid != fold
    fold_tokens <- hal4_subset_tokens(tokens, training)
    fit <- hal4_fit_tokens(
      fold_tokens, y[training], offset[training], lambda = full$lambda,
      penalty_factor = penalty_factor, intercept = intercept, ...)
    validation <- which(!training)
    mu <- vapply(full$lambda, function(value) {
      hal3_predict_tokens(
        fit, data[validation, , drop = FALSE], s = value,
        type = "response", newoffset = offset[validation])
    }, numeric(length(validation)))
    yy <- y[validation]
    cvraw[fold, ] <- colMeans(
      mu - yy * log(pmax(mu, .Machine$double.xmin)))
  }
  full$cvm <- colMeans(cvraw)
  full$cvsd <- apply(cvraw, 2L, stats::sd) / sqrt(nfolds)
  best <- which.min(full$cvm)
  eligible <- which(full$cvm <= full$cvm[best] + full$cvsd[best])
  full$lambda.min <- full$lambda[best]
  full$lambda.1se <- max(full$lambda[eligible])
  full$selected_lambda <- full$lambda.min
  full$foldid <- foldid
  full$cvraw <- cvraw
  class(full) <- c("cv_hal_token_irls_fit", "hal_token_irls_fit")
  full
}
