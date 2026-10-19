#' Fit a sparse Poisson lasso path with custom IRLS
#'
#' Fits an L1-penalized Poisson generalized linear model using outer IRLS and
#' inner coordinate descent directly on a compressed-column sparse matrix.
#' The loss is divided by `nrow(X)`. By default, penalty factors are rescaled to
#' sum to `ncol(X)`, matching `glmnet`.
#'
#' @param X A `Matrix::dgCMatrix`.
#' @param y Nonnegative numeric response.
#' @param offset Numeric observation-level offset.
#' @param lambda Decreasing nonnegative lambda sequence. If `NULL`, a sequence
#'   is generated after a custom-solver pilot fit of unpenalized terms.
#' @param penalty_factor Nonnegative factor for each column.
#' @param intercept Include an unpenalized intercept.
#' @param tol Convergence tolerance.
#' @param maxit_irls Maximum IRLS iterations per lambda.
#' @param maxit_cd Maximum coordinate-descent sweeps per IRLS iteration.
#' @param active_set Use restricted active-set sweeps followed by full KKT checks.
#' @param rescale_penalty_factor Rescale factors to sum to the number of columns,
#'   as `glmnet` does.
#' @param nlambda Number of automatically generated lambdas.
#' @param lambda_min_ratio Smallest/largest lambda ratio for automatic paths.
#' @param check_interrupt_every Frequency of user-interrupt checks.
#'
#' @return A `hal_poisson_irls_fit` object with coefficients, diagnostics, and
#'   component timings for the complete path.
#' @export
hal_poisson_irls <- function(X,
                             y,
                             offset = rep(0, length(y)),
                             lambda = NULL,
                             penalty_factor = rep(1, ncol(X)),
                             intercept = TRUE,
                             tol = 1e-8,
                             maxit_irls = 100L,
                             maxit_cd = 10000L,
                             active_set = TRUE,
                             rescale_penalty_factor = TRUE,
                             nlambda = 100L,
                             lambda_min_ratio = if (nrow(X) < ncol(X)) 0.01 else 1e-4,
                             check_interrupt_every = 25L) {
  if (!inherits(X, "dgCMatrix")) {
    X <- methods::as(X, "dgCMatrix")
  }
  y <- as.numeric(y)
  offset <- as.numeric(offset)
  penalty_factor <- as.numeric(penalty_factor)

  call_solver <- function(lambda_use) {
    hal_poisson_irls_cpp(
      X = X,
      y = y,
      offset = offset,
      lambda = as.numeric(lambda_use),
      penalty_factor = penalty_factor,
      intercept = isTRUE(intercept),
      tol = as.numeric(tol),
      maxit_irls = as.integer(maxit_irls),
      maxit_cd = as.integer(maxit_cd),
      active_set = isTRUE(active_set),
      rescale_penalty_factor = isTRUE(rescale_penalty_factor),
      check_interrupt_every = as.integer(check_interrupt_every)
    )
  }

  if (is.null(lambda) || !length(lambda)) {
    if (length(nlambda) != 1L || is.na(nlambda) || nlambda < 1L ||
        length(lambda_min_ratio) != 1L || !is.finite(lambda_min_ratio) ||
        lambda_min_ratio <= 0 || lambda_min_ratio > 1) {
      stop("Invalid automatic lambda-sequence controls.", call. = FALSE)
    }
    # A very large penalty fixes penalized coordinates at zero while the same
    # IRLS implementation optimizes the intercept and zero-penalty columns.
    pilot <- call_solver(1e30)
    eta <- as.numeric(offset + pilot$a0[1L] + X %*% pilot$beta[, 1L])
    mu <- exp(eta)
    scale <- if (isTRUE(rescale_penalty_factor)) {
      ncol(X) / sum(penalty_factor)
    } else {
      1
    }
    pf_effective <- penalty_factor * scale
    score <- abs(as.numeric(Matrix::crossprod(X, y - mu))) / nrow(X)
    eligible <- pf_effective > 0
    lambda_max <- if (any(eligible)) max(score[eligible] / pf_effective[eligible]) else 0
    if (!is.finite(lambda_max) || lambda_max <= 0) lambda_max <- 1e-8
    lambda <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
                      length.out = as.integer(nlambda)))
  }

  fit <- call_solver(lambda)
  rownames(fit$beta) <- colnames(X)
  colnames(fit$beta) <- format(fit$lambda, digits = 10L)
  fit$call <- match.call()
  fit$dim <- dim(X)
  fit$nnz <- length(X@x)
  fit$offset_supplied <- !missing(offset)
  class(fit) <- "hal_poisson_irls_fit"
  fit
}

.hal_irls_lambda_index <- function(object, s = NULL) {
  if (is.null(s)) {
    if (!is.null(object$selected_lambda)) s <- object$selected_lambda
    else return(seq_along(object$lambda))
  }
  if (length(s) != 1L || !is.finite(s)) stop("`s` must be one finite lambda.")
  which.min(abs(object$lambda - as.numeric(s)))
}

#' @export
coef.hal_poisson_irls_fit <- function(object, s = NULL, ...) {
  index <- .hal_irls_lambda_index(object, s)
  out <- rbind(`(Intercept)` = object$a0[index], object$beta[, index, drop = FALSE])
  colnames(out) <- format(object$lambda[index], digits = 10L)
  out
}

#' @export
predict.hal_poisson_irls_fit <- function(object,
                                         newx,
                                         s = NULL,
                                         type = c("link", "response"),
                                         newoffset = NULL,
                                         ...) {
  type <- match.arg(type)
  if (!inherits(newx, "Matrix")) newx <- methods::as(as.matrix(newx), "dgCMatrix")
  index <- .hal_irls_lambda_index(object, s)
  eta <- sweep(as.matrix(newx %*% object$beta[, index, drop = FALSE]),
               2L, object$a0[index], `+`)
  if (!is.null(newoffset)) {
    if (length(newoffset) != nrow(newx)) stop("`newoffset` has the wrong length.")
    eta <- eta + as.numeric(newoffset)
  }
  if (type == "response") exp(eta) else eta
}

#' @export
coef.cv_hal_poisson_irls_fit <- function(object, s = "lambda.min", ...) {
  if (is.character(s)) s <- object[[s]]
  coef.hal_poisson_irls_fit(object, s = s, ...)
}

#' @export
predict.cv_hal_poisson_irls_fit <- function(object,
                                            newx,
                                            s = "lambda.min",
                                            type = c("link", "response"),
                                            newoffset = NULL,
                                            ...) {
  if (is.character(s)) s <- object[[s]]
  predict.hal_poisson_irls_fit(object, newx = newx, s = s, type = type,
                               newoffset = newoffset, ...)
}

cv_hal_poisson_irls <- function(X, y, offset, lambda = NULL,
                                penalty_factor = rep(1, ncol(X)), nfolds = 5L,
                                foldid = NULL, ...) {
  n <- nrow(X)
  if (is.null(foldid)) foldid <- sample(rep(seq_len(nfolds), length.out = n))
  foldid <- as.integer(foldid)
  if (length(foldid) != n || !setequal(unique(foldid), seq_len(nfolds)))
    stop("`foldid` must assign every row to one of all requested folds.")
  full <- hal_poisson_irls(X, y, offset, lambda = lambda,
                           penalty_factor = penalty_factor, ...)
  cvraw <- matrix(NA_real_, nfolds, length(full$lambda))
  for (fold in seq_len(nfolds)) {
    training <- foldid != fold
    fit <- hal_poisson_irls(X[training, , drop = FALSE], y[training], offset[training],
                            lambda = full$lambda, penalty_factor = penalty_factor, ...)
    mu <- predict(fit, X[!training, , drop = FALSE], newoffset = offset[!training],
                  type = "response")
    yy <- y[!training]
    contribution <- mu - yy * log(pmax(mu, .Machine$double.xmin))
    cvraw[fold, ] <- colMeans(contribution)
  }
  full$cvm <- colMeans(cvraw)
  full$cvsd <- apply(cvraw, 2L, stats::sd) / sqrt(nfolds)
  best <- which.min(full$cvm)
  eligible <- which(full$cvm <= full$cvm[best] + full$cvsd[best])
  full$lambda.min <- full$lambda[best]
  full$lambda.1se <- max(full$lambda[eligible])
  full$foldid <- foldid
  full$cvraw <- cvraw
  class(full) <- c("cv_hal_poisson_irls_fit", "hal_poisson_irls_fit")
  full
}
