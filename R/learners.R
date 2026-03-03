#' Penalized Poisson learner via `glmnet`
#'
#' `Learner_glmnet` is a Reference Class implementing the learner interface
#' used by [Superlearner()] and [fit_learner()].
#'
#' **User-facing API:** users are expected to **initialize** the learner (i.e.,
#' call `Learner_glmnet(...)`) and pass the resulting object to
#' [Superlearner()] or [fit_learner()]. The remaining methods documented below
#' (e.g., `private_fit()`, `private_predictor()`) are part of the internal learner
#' interface and are **not meant to be called directly by users**.
#'
#' **Wrapper role:** this class is a user-friendly wrapper around the existing
#' `glmnet` implementation. The package-specific contribution is to provide a
#' piecewise-constant hazard workflow: create the long-format Poisson data with
#' offsets for time at risk, add the interval ("node") structure for the baseline
#' hazard when requested, and forward standard `glmnet` arguments supplied at
#' initialization to the backend fitter.
#'
#' @section Model:
#' Let \eqn{0=t_0 < t_1 < \cdots < t_m}{0=t0 < t1 < ... < tm} denote time knots and
#' define interval indicators \eqn{I_k(t)=1\{t\in(t_k,t_{k+1}]\}}{Ik(t)=I(t in (tk, t{k+1}])}.
#' The piecewise-constant hazard model is
#' \deqn{
#'   \lambda(t \mid x) = \sum_{k=0}^{m} I_k(t)\,\lambda_k(x),
#'   \qquad \lambda_k(x) = \exp(\beta^\top x + \gamma_k).
#' }{
#'   lambda(t|x) = sum_{k=0}^m Ik(t) * lambda_k(x),
#'   lambda_k(x) = exp(beta^T x + gamma_k).
#' }
#' Penalization is applied to the regression coefficients through the `glmnet`
#' elastic-net penalty. If you want node (baseline) terms to be unpenalized, use
#' `penalty.factor` via `...` (and set it consistently with how your design matrix
#' encodes nodes).
#'
#' @section Fields:
#' \describe{
#'   \item{`covariates` (`character`)}{Names of covariate columns used in the model.}
#'   \item{`cross_validation` (`logical`)}{If `TRUE`, chooses `lambda` by `glmnet::cv.glmnet`.}
#'   \item{`intercept` (`logical`)}{Whether to include an intercept in the backend fit.}
#'   \item{`add_nodes` (`logical`)}{If `TRUE`, include interval ("node") effects encoding the baseline hazard.}
#'   \item{`lambda` (`numeric`)}{If `cross_validation=FALSE`, the `lambda` used in the final fit.}
#'   \item{`formula` (`character`)}{Formula string used to create the design matrix in long format.}
#'   \item{`learner` (`function`)}{Backend fitter (`glmnet::glmnet` or `glmnet::cv.glmnet`).}
#'   \item{`fit_arguments` (`list`)}{Additional arguments forwarded to the backend fitter.}
#' }
#'
#' @section Methods (internal learner interface):
#' \describe{
#'   \item{`initialize(...)`}{Construct and configure the learner. This is the only method users should call.}
#'   \item{`private_fit(data, ...)`}{Internal. Fits a Poisson model with offset `log(tij)` on long-format data.}
#'   \item{`private_predictor(model, newdata, ...)`}{Internal. Predicts hazards on the response scale for long-format `newdata`.}
#' }
#'
#' @examples
#' lrn <- Learner_glmnet(covariates = c("age", "sex"), alpha = 1, cross_validation = TRUE)
#'
#' @export Learner_glmnet
#' @exportClass Learner_glmnet
Learner_glmnet <- setRefClass(
  "Learner_glmnet",
  fields = list(
    covariates = "character",
    cross_validation = "logical",
    intercept = "logical",
    formula = "character",
    learner = "function",
    add_nodes = "logical",
    lambda = "numeric",
    fit_arguments = "list"
  ),
  methods = list(
    initialize = function(covariates = NA_character_,
                          cross_validation = FALSE,
                          intercept = TRUE,
                          add_nodes = TRUE,
                          lambda = NA_real_,
                          ...) {
      .self$covariates <- covariates

      .self$cross_validation <- cross_validation

      .self$intercept <- intercept

      .self$add_nodes <- add_nodes

      tmp <- lambda
      if (is.null(lambda))
        tmp <- NA_real_
      if (length(lambda) == 0L)
        tmp <- NA_real_
      tmp <- as.numeric(tmp)

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula_glmnet(
        covariates = .self$covariates,
        add_nodes = .self$add_nodes
      )



      .self$fit_arguments <- list(...)

      .self$fit_arguments[['family']] <- "poisson"

      .self$fit_arguments[['intercept']] <- .self$intercept

      if (.self$cross_validation) {
        .self$learner = glmnet::cv.glmnet

      } else{
        .self$learner = glmnet::glmnet
        .self$fit_arguments[['lambda']] <- tmp
      }


    },

    private_fit = function(data, ...) {
      .extract_symbols <- function(term) {
        # try as a bare expression, then as a RHS of a formula
        out <- tryCatch(
          all.vars(str2lang(term)),
          error = function(e) {
            tryCatch(
              all.vars(stats::terms(stats::as.formula(
                paste("~", term)
              ))),
              error = function(e2)
                character(0)
            )
          }
        )
        out
      }


      group_cols <- .self$covariates[complete.cases(.self$covariates)]


      group_cols <- group_cols[is.character(group_cols) &
                                 !is.na(group_cols) & nzchar(group_cols)]
      group_cols <- unlist(lapply(group_cols, .extract_symbols), use.names = FALSE)

      group_cols <- unique(group_cols)

      data <- data[, .(tij = sum(tij), deltaij = sum(deltaij)), by = c(group_cols, "node", "k")]

      data <- data[complete.cases(data), ]

      x = Matrix::sparse.model.matrix(formula(.self$formula), data, contrasts.arg = NULL)[, -1]




      pf <- 1 - (grepl("node", colnames(x)))


      if (.self$cross_validation) {
        # Cross-validation is done using cv.glmnet's internally generated
        # lambda path. We then refit a glmnet model at the lambda that minimizes
        # the (Poisson) deviance, matching the approach used for Learner_hal.

        cv_args <- .self$fit_arguments

        # Ensure deviance is used unless the user explicitly provided another
        # type.measure.
        if (is.null(cv_args[['type.measure']])) {
          cv_args[['type.measure']] <- 'deviance'
        }

        cv_fit <- tryCatch(
          suppressWarnings(
            do.call(glmnet::cv.glmnet, c(
              list(
                x = x,
                y = as.vector(data[["deltaij"]]),
                offset = log(as.vector(data[["tij"]])),
                penalty.factor =pf
              ),
              cv_args
            ))
          ),
          error = function(e) {
            suppressWarnings(
              do.call(glmnet::cv.glmnet, c(
                list(
                  x = x,
                  y = as.vector(data[["deltaij"]]),
                  offset = log(as.vector(data[["tij"]])),
                  penalty.factor = rep(1, ncol(x))
                ),
                cv_args
              ))
            )
          }
        )

        # suppressWarnings(cv_fit <- do.call(cv.glmnet, c(
        #   list(
        #     x = x,
        #     y = as.vector(data[["deltaij"]]),
        #     offset = log(as.vector(data[["tij"]]))
        #   ),cv_args
        # ))


        # lambda.min corresponds to the minimizer of cvm (expected deviance)
        # over cv.glmnet's lambda path.
        lambda_opt <- cv_fit$lambda.min

        # Store selected lambda for potential downstream use.
        .self$lambda <- as.numeric(lambda_opt)

        # Refit glmnet at the selected lambda.
        glmnet_args <- .self$fit_arguments
        glmnet_args[['lambda']] <- as.numeric(lambda_opt)

        # Remove CV-only arguments if present.
        glmnet_args[['nfolds']] <- NULL
        glmnet_args[['foldid']] <- NULL
        glmnet_args[['type.measure']] <- NULL
        glmnet_args[['keep']] <- NULL
        glmnet_args[['grouped']] <- NULL
        glmnet_args[['parallel']] <- NULL
        glmnet_args[['penalty.factor']] <- pf


        suppressWarnings(out <- do.call(glmnet::glmnet, c(glmnet_args, list(
          x = x,
          y = as.numeric(data[["deltaij"]]),
          offset = log(data[["tij"]])
        ))))

      }else {
        suppressWarnings(out <- do.call(.self$learner, c(
          .self$fit_arguments,
          list(
            x = x,
            y = as.numeric(data[["deltaij"]]),
            offset = log(data[["tij"]])
          )
        )))
      }

      return(out)

    },

    private_predictor = function(model, newdata, ...) {
      is_empty_model <- function(model) {
        if (is.null(model))
          return(TRUE)
        if (inherits(model, "cv.glmnet")) {
          if (is.null(model$glmnet.fit))
            return(TRUE)
          beta <- model$glmnet.fit$beta
        } else {
          beta <- model$beta
        }
        if (is.null(beta))
          return(TRUE)
        if (!is.matrix(beta) &&
            !inherits(beta, "Matrix"))
          return(FALSE)
        nrow(beta) == 0 || ncol(beta) == 0
      }

      if (is_empty_model(model)) {
        return(rep(NA_real_, nrow(newdata)))
      }

      out <- predict(
        model,
        newx = Matrix::sparse.model.matrix(formula(.self$formula), newdata, contrasts.arg = NULL)[, -1],
        newoffset = log(1),
        type = "response",
        ...
      )

      return(out)


    }

  )
)

#' HAL learner for piecewise Poisson hazards
#'
#' `Learner_hal` is a Reference Class implementing the learner interface
#' used by [Superlearner()] and [fit_learner()].
#'
#' **User-facing API:** users should **only initialize** the learner and pass it
#' to [Superlearner()] / [fit_learner()]. The methods `private_fit()` and
#' `private_predictor()` (and any basis-construction helpers) are part of the
#' internal learner interface and are **not meant to be called directly by users**.
#'
#' **Wrapper role:** this class provides a piecewise-constant hazard wrapper around
#' a HAL-style indicator-basis construction, estimated by L1-penalized Poisson
#' regression using a `glmnet` backend. The package-specific contribution is to
#' (i) construct the long-format Poisson representation with offsets for time at
#' risk, (ii) generate indicator bases compatible with piecewise hazards, and
#' (iii) forward backend fitting arguments supplied via `...`.
#'
#' @section Model:
#' Let \eqn{0=t_0 < t_1 < \cdots < t_m}{0=t0 < t1 < ... < tm} denote time knots and
#' define interval indicators \eqn{I_k(t)=1\{t\in(t_k,t_{k+1}]\}}{Ik(t)=I(t in (tk, t{k+1}])}.
#' The HAL piecewise-constant hazard model is
#' \deqn{
#'   \lambda(t \mid x) = \sum_{k=0}^{m} I_k(t)\,\exp\{f(t,x)\},
#' }{
#'   lambda(t|x) = sum_{k=0}^m Ik(t) * exp(f(t,x)).
#' }
#' where \eqn{f(t,x)}{f(t,x)} is approximated by a finite linear combination of
#' indicator basis functions.
#'
#' @section Two-covariate illustration:
#' Let \eqn{x=(x_1,x_2)}{x=(x1,x2)} be two covariates and let
#' \eqn{t_0 < t_1 < \cdots < t_R}{t0 < t1 < ... < tR} be time grid points used to
#' create step functions in time. Choose covariate cutpoints
#' \eqn{c_{1,1},\ldots,c_{1,K_1}}{c1_1,...,c1_K1} for \eqn{x_1}{x1} and
#' \eqn{c_{2,1},\ldots,c_{2,K_2}}{c2_1,...,c2_K2} for \eqn{x_2}{x2}.
#'
#' Define indicator bases:
#' \deqn{B_r(t) = 1\{t_r \le t\}}{Br(t) = I(tr <= t),}
#' \deqn{B_{1,p}(x) = 1\{c_{1,p} \le x_1\}}{B1p(x) = I(c1p <= x1),}
#' \deqn{B_{2,q}(x) = 1\{c_{2,q} \le x_2\}}{B2q(x) = I(c2q <= x2).}
#'
#' A main-effects HAL approximation on the log-hazard scale can be written as:
#' \deqn{
#'   f_\beta(t,x) = \beta_0
#'   + \sum_{r=1}^R \beta_r B_r(t)
#'   + \sum_{r=1}^R\sum_{p=1}^{K_1} \beta_{r,1,p} B_r(t) B_{1,p}(x)
#'   + \sum_{r=1}^R\sum_{q=1}^{K_2} \beta_{r,2,q} B_r(t) B_{2,q}(x).
#' }{
#'   f_beta(t,x) = beta0
#'   + sum_{r=1}^R beta_r Br(t)
#'   + sum_{r=1}^R sum_{p=1}^K1 beta_{r,1,p} Br(t) B1p(x)
#'   + sum_{r=1}^R sum_{q=1}^K2 beta_{r,2,q} Br(t) B2q(x).
#' }
#' If `max_degree >= 2`, the learner additionally includes interaction bases such as
#' \deqn{
#'   \sum_{r=1}^R\sum_{p=1}^{K_1}\sum_{q=1}^{K_2}
#'   \beta_{r,12,pq} B_r(t) B_{1,p}(x) B_{2,q}(x).
#' }{
#'   sum_{r=1}^R sum_{p=1}^K1 sum_{q=1}^K2 beta_{r,12,pq} Br(t) B1p(x) B2q(x).
#' }
#'
#' @section How reference class parameters map to the model:
#' \describe{
#'   \item{`covariates`}{Covariate columns used to build covariate indicator bases.}
#'   \item{`num_knots`}{Controls the number of cutpoints per covariate used for indicator bases.}
#'   \item{`max_degree`}{Maximum interaction order included in the basis expansion.}
#'   \item{`add_nodes`}{If `TRUE`, includes interval ("node") structure for the baseline hazard.}
#'   \item{`intercept`}{Whether the backend penalized regression includes an intercept term.}
#'   \item{`cross_validation`}{If `TRUE`, selects the penalty level using `glmnet::cv.glmnet`.}
#'   \item{`fit_arguments`}{Additional arguments forwarded to the `glmnet` backend (e.g. `nfolds`).}
#' }
#'
#' @section Fields:
#' \describe{
#'   \item{`covariates` (`character`)}{Names of covariate columns used in the basis.}
#'   \item{`cross_validation` (`logical`)}{Whether to use `cv.glmnet` to select the penalty.}
#'   \item{`intercept` (`logical`)}{Backend intercept flag.}
#'   \item{`add_nodes` (`logical`)}{Whether node (time-interval) effects are included.}
#'   \item{`max_degree` (`integer`)}{Maximum interaction order.}
#'   \item{`num_knots` (`numeric`)}{Knots used for basis construction.}
#'   \item{`lambda_opt` (`numeric`)}{Selected penalty level when using cross-validation.}
#'   \item{`fit_arguments` (`list`)}{Extra backend arguments forwarded to `glmnet`.}
#' }
#'
#' @section Methods (internal learner interface):
#' \describe{
#'   \item{`initialize(...)`}{Construct and configure the learner. This is the only method users should call.}
#'   \item{`private_fit(data, ...)`}{Internal. Builds bases and fits the penalized Poisson model with offset `log(tij)`.}
#'   \item{`private_predictor(model, newdata, ...)`}{Internal. Evaluates the fitted approximation and returns hazards on the response scale.}
#' }
#'
#' @examples
#' lrn <- Learner_hal(covariates = c("age", "sex"), max_degree = 2L, num_knots = c(10L, 5L))
#'
#' @export Learner_hal
#' @exportClass Learner_hal
Learner_hal <- setRefClass(
  "Learner_hal",
  fields = list(
    covariates = "character",
    cross_validation = "logical",
    intercept = "logical",
    formula = "character",
    learner = "function",
    add_nodes = "logical",
    max_degree = "integer",
    lambda_opt = "numeric",
    maxit_prefit = "numeric",
    fit_arguments = "list",
    covariates_attributes_matrix = "list",
    num_knots = "numeric"
  ),
  methods = list(
    initialize = function(covariates = NA_character_,
                          intercept = FALSE,
                          cross_validation = TRUE,
                          add_nodes = TRUE,
                          max_degree = 2,
                          maxit_prefit = NA_real_,
                          num_knots = c(10L, 5L),
                          ...) {
      .self$covariates <- covariates

      .self$cross_validation <- cross_validation

      .self$intercept <- intercept

      .self$add_nodes <- add_nodes

      .self$max_degree <- max_degree

      .self$num_knots <- num_knots

      # normalize user input
      # tmp <- lambda_grid
      # if (is.null(lambda_grid))
      #   tmp <- NA_real_
      # if (length(lambda_grid) == 0L)
      #   tmp <- NA_real_
      #
      # .self$lambda_grid <- tmp

      .self$maxit_prefit <- maxit_prefit

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula_hal(
        covariates = .self$covariates,
        intercept = FALSE,
        #in the glmnet case, intercept is handled separately.
        add_nodes = .self$add_nodes
      )

      if (.self$cross_validation) {
        .self$learner <- glmnet::cv.glmnet
      } else {
        .self$learner <- glmnet::glmnet
      }

      .self$fit_arguments <- list(...)

      .self$covariates_attributes_matrix <- list()

      .self$fit_arguments[['family']] <- "poisson"

      .self$fit_arguments[['intercept']] <- .self$intercept




    },
    hal_basis = function(vars,
                         DT,
                         max_interaction = 2L,
                         knots_per_order = rep(10L, 2L)) {
      stopifnot(requireNamespace("data.table", quietly = TRUE))
      stopifnot(requireNamespace("Matrix", quietly = TRUE))

      if (!data.table::is.data.table(DT))
        DT <- data.table::as.data.table(DT)
      if (length(knots_per_order) != max_interaction)
        stop("knots_per_order must have length max_interaction")
      if (!all(vars %in% names(DT)))
        stop("vars not all found in DT")

      n <- nrow(DT)

      DT_local <- DT[, ..vars]
      for (v in vars) {
        x <- DT_local[[v]]
        if (is.character(x) ||
            is.logical(x))
          DT_local[[v]] <- factor(x)
      }

      primitives <- setNames(vector("list", length(vars)), vars)
      meta_per_var <- setNames(vector("list", length(vars)), vars)

      for (v in vars) {
        res <- mk_main(DT_local[[v]], v, K = knots_per_order[1L])
        primitives[[v]] <- res
        meta_per_var[[v]] <- res$var_meta
      }

      # state for chunked construction (purely threaded)
      state <- list(
        I_chunks = vector("list", 0L),
        J_chunks = vector("list", 0L),
        name_chunks = vector("list", 0L),
        colmeta_chunks = vector("list", 0L),
        chunk_used = 0L,
        p = 0L
      )

      # mains
      for (v in vars) {
        prim <- primitives[[v]]
        if (!length(prim$idxs))
          next
        mains_meta <- lapply(prim$prim_meta, function(m)
          c(list(type = "main"), m))
        state <- add_cols(state, prim$idxs, prim$names, mains_meta)
      }

      # interactions
      if (max_interaction >= 2L) {
        for (d in 2:max_interaction) {
          if (length(vars) < d)
            break

          cmbs <- utils::combn(vars, d, simplify = FALSE)
          Kcap <- knots_per_order[d]

          for (cmb in cmbs) {
            prim_d <- lapply(cmb, function(v) {
              idxs <- primitives[[v]]$idxs
              pm   <- primitives[[v]]$prim_meta
              if (!length(idxs))
                return(list(
                  idxs = list(),
                  names = character(),
                  meta = list()
                ))

              is_num <- length(pm) > 0L &&
                identical(pm[[1]]$kind, "numeric")
              ncap <- if (is_num)
                min(length(idxs), Kcap)
              else
                length(idxs)

              list(idxs  = idxs[seq_len(ncap)],
                   names = primitives[[v]]$names[seq_len(ncap)],
                   meta  = pm[seq_len(ncap)])
            })

            if (any(vapply(prim_d, function(z)
              length(z$idxs), integer(1L)) == 0L))
              next

            lens <- vapply(prim_d, function(z)
              length(z$idxs), integer(1L))
            total <- prod(lens)
            if (!total)
              next

            for (t in 0:(total - 1L)) {
              sel <- integer(d)
              rem <- t
              for (j in seq_len(d)) {
                sel[j] <- (rem %% lens[j]) + 1L
                rem <- rem %/% lens[j]
              }

              chosen_idx <- vector("list", d)
              nm_parts <- character(d)
              parts_meta <- vector("list", d)

              for (j in seq_len(d)) {
                chosen_idx[[j]] <- prim_d[[j]]$idxs[[sel[j]]]
                nm_parts[j] <- prim_d[[j]]$names[[sel[j]]]
                parts_meta[[j]] <- prim_d[[j]]$meta[[sel[j]]]
              }

              idx_prod <- interN_cpp(chosen_idx)
              if (!length(idx_prod))
                next

              state <- add_cols(
                state,
                idxs_list = list(idx_prod),
                nm_vec = paste0(nm_parts, collapse = ":"),
                col_meta_list = list(
                  list(
                    type = "interaction",
                    order = d,
                    parts = parts_meta
                  )
                )
              )
            }
          }
        }
      }

      if (state$chunk_used == 0L)
        stop("No basis columns created.")

      I_chunks <- state$I_chunks
      J_chunks <- state$J_chunks
      name_chunks <- state$name_chunks
      colmeta_chunks <- state$colmeta_chunks

      i_vec <- unlist(I_chunks, use.names = FALSE)
      j_vec <- unlist(J_chunks, use.names = FALSE)

      names_all <- unlist(name_chunks, use.names = FALSE)
      per_col_meta <- unlist(colmeta_chunks, recursive = FALSE, use.names = FALSE)

      X <- Matrix::sparseMatrix(
        i = if (length(i_vec))
          i_vec
        else
          integer(),
        j = if (length(j_vec))
          j_vec
        else
          integer(),
        x = if (length(i_vec))
          1L
        else
          integer(),
        dims = c(n, if (length(j_vec))
          max(j_vec)
          else
            0L),
        dimnames = list(NULL, names_all)
      )

      list(
        X = X,
        colnames = names_all,
        meta = list(per_var = meta_per_var, per_col = per_col_meta)
      )
    },

    hal_prepare_new = function(DT_new, hal_obj) {
      stopifnot(requireNamespace("data.table", quietly = TRUE))
      stopifnot(requireNamespace("Matrix", quietly = TRUE))

      if (!data.table::is.data.table(DT_new))
        DT_new <- data.table::as.data.table(DT_new)

      vars <- names(hal_obj$meta$per_var)
      if (!all(vars %in% names(DT_new)))
        stop("DT_new is missing required variables")

      DT_local <- DT_new[, ..vars]

      if (!exists("interN_cpp", mode = "function", inherits = TRUE)) {
        stop("interN_cpp() not found. Ensure the C++ code is compiled and loaded.")
      }

      n <- nrow(DT_local)
      colnames_tr <- hal_obj$colnames
      colmeta <- hal_obj$meta$per_col
      if (length(colmeta) != length(colnames_tr))
        stop("hal_obj meta mismatch")

      # Stable key for numeric cutpoints (must be used consistently)
      key_num <- function(z)
        formatC(z, digits = 17, format = "fg")

      # --- coerce DT_local types using training meta_per_var ---
      for (v in vars) {
        per <- hal_obj$meta$per_var[[v]]

        if (!is.null(per$levels)) {
          # factor variables: enforce training levels
          x <- DT_local[[v]]

          if (is.character(x) || is.logical(x))
            x <- as.character(x)
          if (!is.factor(x))
            x <- as.character(x)

          lv_tr <- per$levels
          DT_local[[v]] <- factor(x, levels = lv_tr)

          if (anyNA(DT_local[[v]])) {
            bad <- unique(as.character(x[is.na(DT_local[[v]])]))
            stop(sprintf(
              "DT_new has unseen level(s) in %s: %s",
              v,
              paste(bad, collapse = ", ")
            ))
          }

        } else if (!is.null(per$cutpoints)) {
          DT_local[[v]] <- as.numeric(DT_local[[v]])
        } else {
          x <- DT_local[[v]]
          if (is.character(x) ||
              is.logical(x))
            DT_local[[v]] <- factor(x)
        }
      }
      # -------------------------------------------------------------

      # Precompute which primitives are needed (from per_col)
      need_num <- lapply(vars, function(v)
        numeric())
      names(need_num) <- vars
      need_fac <- lapply(vars, function(v)
        character())
      names(need_fac) <- vars

      for (m in colmeta) {
        if (m$type == "main") {
          if (m$kind == "numeric")
            need_num[[m$var]] <- unique(c(need_num[[m$var]], m$cutpoint))
          else
            need_fac[[m$var]] <- unique(c(need_fac[[m$var]], m$level))
        } else if (m$type == "interaction") {
          for (p in m$parts) {
            if (p$kind == "numeric")
              need_num[[p$var]] <- unique(c(need_num[[p$var]], p$cutpoint))
            else
              need_fac[[p$var]] <- unique(c(need_fac[[p$var]], p$level))
          }
        } else {
          stop("Unknown column meta type")
        }
      }

      # Build index maps: var -> key -> sorted row indices
      num_map <- vector("list", length(vars))
      names(num_map) <- vars
      fac_map <- vector("list", length(vars))
      names(fac_map) <- vars

      for (v in vars) {
        x <- DT_local[[v]]

        # numeric: STANDARD HAL uses I(x >= cut)
        if (is.numeric(x) && length(need_num[[v]])) {
          cuts <- sort(unique(need_num[[v]]))
          nn <- which(!is.na(x) & is.finite(x))

          if (length(nn)) {
            ord <- nn[order(x[nn])]
            xs <- x[ord]
            N <- length(ord)

            # k_lt[j] = number of xs < cuts[j]
            # (left.open=TRUE makes equality go to previous interval => strict <)
            k_lt <- findInterval(cuts, xs, left.open = TRUE)

            idxs <- lapply(k_lt, function(k) {
              start <- k + 1L
              if (start <= N)
                sort.int(ord[start:N], method = "quick")
              else
                integer()
            })

          } else {
            idxs <- lapply(cuts, function(.)
              integer())
          }

          names(idxs) <- key_num(cuts)
          num_map[[v]] <- idxs
        }

        # factor
        if (is.factor(x) && length(need_fac[[v]])) {
          chr <- as.character(x)
          lv_needed <- need_fac[[v]]
          idxs <- lapply(lv_needed, function(L)
            which(!is.na(chr) & chr == L))
          names(idxs) <- lv_needed
          fac_map[[v]] <- idxs
        }
      }

      # Assemble sparse triplets in training column order
      I_chunks <- vector("list", length(colmeta))
      J_chunks <- vector("list", length(colmeta))

      for (j in seq_along(colmeta)) {
        m <- colmeta[[j]]

        if (m$type == "main") {
          if (m$kind == "numeric") {
            idx <- num_map[[m$var]][[key_num(m$cutpoint)]]
          } else {
            idx <- fac_map[[m$var]][[m$level]]
          }

        } else {
          # interaction
          parts_idx <- vector("list", length(m$parts))
          ok <- TRUE

          for (k in seq_along(m$parts)) {
            p <- m$parts[[k]]
            if (p$kind == "numeric") {
              parts_idx[[k]] <- num_map[[p$var]][[key_num(p$cutpoint)]]
            } else {
              parts_idx[[k]] <- fac_map[[p$var]][[p$level]]
            }
            if (is.null(parts_idx[[k]]) ||
                !length(parts_idx[[k]])) {
              ok <- FALSE
              break
            }
          }

          if (!ok) {
            idx <- integer()
          } else if (length(parts_idx) == 1L) {
            idx <- parts_idx[[1L]]
          } else {
            idx <- interN_cpp(parts_idx)
          }
        }

        if (length(idx)) {
          I_chunks[[j]] <- idx
          J_chunks[[j]] <- rep.int(j, length(idx))
        } else {
          I_chunks[[j]] <- integer()
          J_chunks[[j]] <- integer()
        }
      }

      i_vec <- unlist(I_chunks, use.names = FALSE)
      j_vec <- unlist(J_chunks, use.names = FALSE)

      Matrix::sparseMatrix(
        i = if (length(i_vec))
          i_vec
        else
          integer(),
        j = if (length(j_vec))
          j_vec
        else
          integer(),
        x = if (length(i_vec))
          1L
        else
          integer(),
        dims = c(n, length(colnames_tr)),
        dimnames = list(NULL, colnames_tr)
      )
    },



    private_fit = function(data, ...) {
      data_copy = data.table::copy(data)

      group_cols <- .self$covariates[complete.cases(.self$covariates)]

      data_copy <- data_copy[, .(tij = sum(tij), deltaij = sum(deltaij)), by = c(group_cols, "node", "k")]

      data_copy <- data_copy[complete.cases(data_copy), ]

      x_pp <- hal_basis(
        vars = c(group_cols, "node"),
        DT = data_copy,
        max_interaction = .self$max_degree,
        knots_per_order = .self$num_knots
      )

      pf <- 1 - grepl('^I\\(\\s*node\\s*==[^)]*\\)$', x_pp$colnames)


      .self$covariates_attributes_matrix <- x_pp[c("colnames", "meta")]



      if (.self$cross_validation) {
        ## 1) Prefit with cv.glmnet to obtain the lambda path + CV curve
        prefit_args <- .self$fit_arguments

        if (!is.na(.self$maxit_prefit)) {
          prefit_args[["maxit"]] <- .self$maxit_prefit
        }

        cv_fit <- tryCatch(
          suppressWarnings(
            do.call(glmnet::cv.glmnet, c(
              list(
                x      = x_pp$X,
                y      = as.numeric(data_copy[["deltaij"]]),
                offset = log(data_copy[["tij"]]),
                penalty.factor =pf
              ),
              prefit_args
            ))
          ),
          error = function(e) {
            suppressWarnings(
              do.call(glmnet::cv.glmnet, c(
                list(
                  x = x,
                  y = as.vector(data[["deltaij"]]),
                  offset = log(as.vector(data[["tij"]])),
                  penalty.factor = rep(1, ncol(x))
                ),
                prefit_args
              ))
            )
          }
        )





        ## 2) Choose lambda that minimizes cvm (this is cv_fit$lambda.min)
        lambda_min_cvm <- cv_fit$lambda.min
        .self$lambda_opt <- lambda_min_cvm

        ## 3) Refit a plain glmnet at that single lambda (no CV)
        glmnet_args <- .self$fit_arguments
        glmnet_args[["lambda"]] <- lambda_min_cvm
        glmnet_args[["nfolds"]] <- NULL  # not used by glmnet, but remove if present
        glmnet_args[["penalty.factor"]] <- pf


        suppressWarnings(fit <- do.call(glmnet::glmnet, c(
          glmnet_args, list(
            x      = x_pp$X,
            y      = as.numeric(data_copy[["deltaij"]]),
            offset = log(data_copy[["tij"]])
          )
        )))

       } else {
        glmnet_args <- .self$fit_arguments

        # If user supplied lambda_grid (single lambda or vector), use it directly
        # if (!all(is.na(.self$lambda_grid))) {
        #   lambda_user <- as.numeric(.self$lambda_grid)
        #   lambda_user <- lambda_user[complete.cases(lambda_user)]
        #
        #   if (length(lambda_user) > 0L) {
        #     glmnet_args[["lambda"]] <- lambda_user
        #
        #     # store lambda_opt only if a single lambda was provided
        #     if (length(lambda_user) == 1L) {
        #       .self$lambda_opt <- lambda_user
        #     }
        #   }
        # }

        suppressWarnings(fit <- do.call(glmnet::glmnet, c(
          glmnet_args, list(
            x      = x_pp$X,
            y      = as.numeric(data_copy[["deltaij"]]),
            offset = log(data_copy[["tij"]])
          )
        )))
      }

      return(fit)

    },

    private_predictor = function(model, newdata, ...) {
      is_empty_model <- function(fit) {
        if (is.null(fit))
          return(TRUE)
        df <- if (inherits(fit, "cv.glmnet"))
          fit$glmnet.fit$df
        else
          fit$df
        if (!is.null(df))
          return(all(df == 0))
        beta <- if (inherits(fit, "cv.glmnet"))
          fit$glmnet.fit$beta
        else
          fit$beta
        if (is.null(beta))
          return(TRUE)
        all(beta == 0)
      }

      if (is_empty_model(model)) {
        return(rep(NA_real_, nrow(newdata)))
      }

      X_new <- .self$hal_prepare_new(newdata, .self$covariates_attributes_matrix)

      # if(.self$cross_validation){
      out <- predict(
        model,
        newx = X_new,
        newoffset = log(1),
        type = "response",
        ...
      )


      # }

      return(out)


    }



  )
)


#' GAM learner via `mgcv::bam`
#'
#' `Learner_gam` is a Reference Class implementing the learner interface
#' used by [Superlearner()] and [fit_learner()].
#'
#' **User-facing API:** users should **only initialize** the learner and pass it
#' to [Superlearner()] / [fit_learner()]. The methods `private_fit()` and
#' `private_predictor()` are part of the internal learner interface and are
#' **not meant to be called directly by users**.
#'
#' **Wrapper role:** this class wraps `mgcv::bam` in a piecewise-constant hazard
#' workflow. The package-specific contribution is to provide a convenient
#' interface for the long-format Poisson likelihood with offsets for time at risk,
#' and optional node terms encoding the baseline hazard, while forwarding standard
#' `mgcv::bam` arguments supplied via `...`.
#'
#' @section Model:
#' Let \eqn{0=t_0 < t_1 < \cdots < t_m}{0=t0 < t1 < ... < tm} denote time knots and
#' define interval indicators \eqn{I_k(t)=1\{t\in(t_k,t_{k+1}]\}}{Ik(t)=I(t in (tk, t{k+1}])}.
#' The piecewise-constant hazard model with an additive predictor is
#' \deqn{
#'   \lambda(t \mid x) = \sum_{k=0}^{m} I_k(t)\,\exp\{\eta(x) + \gamma_k\}.
#' }{
#'   lambda(t|x) = sum_{k=0}^m Ik(t) * exp(eta(x) + gamma_k).
#' }
#' The additive predictor \eqn{\eta(x)}{eta(x)} is constructed from `covariates`
#' (smooth terms such as `s(age)` and/or linear terms) and estimated by `mgcv`.
#'
#' @param covariates `character`. Right-hand-side terms, including `mgcv` smooths
#'   (e.g. `"s(age)"`) and/or linear terms (e.g. `"value_LDL"`).
#' @param cross_validation `logical`. Included for compatibility with the learner
#'   interface; smoothing selection is controlled by `mgcv` and arguments in `...`.
#'
#' @section Fields:
#' \describe{
#'   \item{`covariates` (`character`)}{Terms used to build the additive predictor (may include `s()` terms).}
#'   \item{`cross_validation` (`logical`)}{Workflow flag; see Details.}
#'   \item{`intercept` (`logical`)}{Whether to include an intercept.}
#'   \item{`add_nodes` (`logical`)}{If `TRUE`, include interval ("node") effects encoding the baseline hazard.}
#'   \item{`formula` (`character`)}{Formula string passed to `mgcv::bam`.}
#'   \item{`learner` (`function`)}{Backend fitter (`mgcv::bam`).}
#'   \item{`fit_arguments` (`list`)}{Additional arguments forwarded to `mgcv::bam`.}
#' }
#'
#' @section Methods (internal learner interface):
#' \describe{
#'   \item{`initialize(...)`}{Construct and configure the learner. This is the only method users should call.}
#'   \item{`private_fit(data, ...)`}{Internal. Fits a Poisson GAM with offset `log(tij)` on long-format data.}
#'   \item{`private_predictor(model, newdata, ...)`}{Internal. Predicts hazards on the response scale.}
#' }
#'
#' @examples
#' lrn <- Learner_gam(covariates = c("s(age)", "value_LDL"))
#'
#' @export Learner_gam
#' @exportClass Learner_gam
Learner_gam <- setRefClass(
  "Learner_gam",
  fields = list(
    covariates = "character",
    cross_validation = "logical",
    intercept = "logical",
    formula = "character",
    learner = "function",
    add_nodes = "logical",
    fit_arguments = "list"
  ),
  methods = list(
    initialize = function(covariates = NULL,
                          cross_validation = FALSE,
                          intercept = TRUE,
                          add_nodes = TRUE,
                          ...) {
      .self$covariates <- covariates
      .self$cross_validation <- cross_validation
      .self$intercept <- intercept
      .self$add_nodes <- add_nodes

      .self$formula <- create_formula_gam(
        covariates = .self$covariates,
        intercept = .self$intercept,
        add_nodes = .self$add_nodes
      )

      .self$learner <- mgcv::bam

      .self$fit_arguments <- list(...)
      .self$fit_arguments[['family']] <- poisson()
    },
    private_fit = function(data, ...) {
      .extract_symbols <- function(term) {
        # try as a bare expression, then as a RHS of a formula
        out <- tryCatch(
          all.vars(str2lang(term)),
          error = function(e) {
            tryCatch(
              all.vars(stats::terms(stats::as.formula(
                paste("~", term)
              ))),
              error = function(e2)
                character(0)
            )
          }
        )
        out
      }


      group_cols <- .self$covariates[complete.cases(.self$covariates)]


      group_cols <- group_cols[is.character(group_cols) &
                                 !is.na(group_cols) & nzchar(group_cols)]
      group_cols <- unlist(lapply(group_cols, .extract_symbols), use.names = FALSE)

      group_cols <- unique(group_cols)

      data <- data[, .(tij = sum(tij), deltaij = sum(deltaij)), by = c(group_cols, "node", "k")]

      data <- data[complete.cases(data), ]


      .self$fit_arguments$formula <- as.formula(.self$formula)
      # .self$fit_arguments$data <- data
      # .self$fit_arguments$offset <- log(data[['tij']])

      fit <- do.call(.self$learner, c(.self$fit_arguments, list(
        data = data, offset = log(data[['tij']])
      )))
      return(fit)
    },

    private_predictor = function(model, newdata, ...) {
      pred <- predict(
        model,
        newdata = newdata[node %in% model$xlevels$node, ],
        type = "response",
        offset = log(1),
        #log(newdata[['tij']]),
        ...
      )


      if (all((levels(newdata$node) %in% model$xlevels$node))) {
        return(pred)

      } else{
        #
        newdata[node %in% model$xlevels$node, predictions_model := pred]

        return(as.array(newdata$predictions_model))




      }
    }
  )
)
