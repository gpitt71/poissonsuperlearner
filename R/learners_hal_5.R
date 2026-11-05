#' HAL learner using `hal9001` basis construction
#'
#' `Learner_hal5` retains the Poisson/glmnet fitting and prediction interface of
#' [Learner_hal()], but delegates enumeration and evaluation of the zero-order
#' HAL basis to `hal9001`.
#'
#' @inheritParams Learner_hal
#' @return An initialized `Learner_hal5` Reference Class object.
#' @export Learner_hal5
#' @exportClass Learner_hal5
Learner_hal5 <- methods::setRefClass(
  "Learner_hal5",
  contains = "Learner_hal",
  methods = list(
    private_fit = function(data, cause, grid_nodes, ...) {
      fit <- tryCatch({
        widths <- c(diff(grid_nodes), 1.0)
        needed_cols <- unique(c(.self$basic_covariates, "node", "tij", "deltaij"))
        train_d <- data.table::copy(data[, .SD, .SDcols = needed_cols])

        if (length(.self$basic_covariates) > 0L) {
          train_d[, gid := .GRP, by = c(.self$basic_covariates)]
          group_key <- unique(
            train_d[, c(.self$basic_covariates, "gid"), with = FALSE]
          )
        } else {
          group_key <- data.table::data.table(gid = 1L)
          train_d[, gid := 1L]
        }

        train_d <- train_d[!is.na(gid) & !is.na(node) & !is.na(tij)]
        train_d <- train_d[, .(
          tij = sum(tij),
          deltaij = sum(as.numeric(deltaij == cause)),
          N_terminal = .N
        ), by = .(gid, node)]

        train_d[, node_ix := match(node, grid_nodes)]
        data.table::setorder(train_d, gid, node_ix)

        ans <- expand_terminal_grouped_cpp(
          gid = train_d$gid,
          node_ix = train_d$node_ix,
          tij = train_d$tij,
          deltaij = train_d$deltaij,
          N = train_d$N_terminal,
          widths = widths
        )

        train_d <- data.table::as.data.table(ans)
        train_d[, node := grid_nodes[node_ix]]
        train_d <- group_key[train_d, on = "gid"]
        data.table::setorder(train_d, gid, node_ix)
        train_d[, node := factor(
          node_ix,
          levels = seq_along(grid_nodes),
          labels = paste0("n", seq_along(grid_nodes))
        )]

        basis_vars <- c(.self$basic_covariates, "node")
        encoding <- setNames(vector("list", length(basis_vars)), basis_vars)
        encoded <- list()

        for (variable in basis_vars) {
          value <- train_d[[variable]]
          if (is.character(value) || is.logical(value)) {
            value <- factor(value)
          }

          if (is.factor(value)) {
            variable_levels <- levels(value)
            encoding[[variable]] <- list(type = "factor", levels = variable_levels)
            if (length(variable_levels) > 1L) {
              for (level in variable_levels[-1L]) {
                column_name <- paste0(variable, "_", make.names(level))
                encoded[[column_name]] <- as.numeric(value == level)
              }
            }
          } else {
            encoding[[variable]] <- list(type = "numeric")
            encoded[[variable]] <- as.numeric(value)
          }
        }

        X_raw <- do.call(cbind, encoded)
        if (is.null(dim(X_raw))) {
          X_raw <- matrix(X_raw, ncol = 1L)
        }
        storage.mode(X_raw) <- "double"

        basis_list <- hal9001::enumerate_basis(
          X_raw,
          max_degree = .self$max_degree,
          smoothness_orders = rep(0L, ncol(X_raw)),
          include_zero_order = FALSE,
          include_lower_order = FALSE,
          num_knots = .self$num_knots
        )
        X_hal <- hal9001::make_design_matrix(X_raw, basis_list)

        if (nrow(X_hal) == 0L || ncol(X_hal) == 0L) {
          return(make_failed_fit("Empty HAL basis in Learner_hal5::private_fit"))
        }

        # Remove constant columns. In particular, hal9001 enumerates the
        # minimum cutoff, whose indicator duplicates the glmnet intercept.
        column_sums <- Matrix::colSums(X_hal)
        keep <- column_sums > 0 & column_sums < nrow(X_hal)
        X_hal <- X_hal[, keep, drop = FALSE]
        basis_list <- basis_list[keep]

        if (ncol(X_hal) == 0L) {
          return(make_failed_fit("Only constant HAL bases in Learner_hal5::private_fit"))
        }

        node_columns <- grepl("^node_", colnames(X_raw))
        penalty_factor <- vapply(basis_list, function(basis) {
          is_node_main <- length(basis$cols) == 1L &&
            node_columns[basis$cols] && basis$cutoffs == 1
          as.numeric(!is_node_main)
        }, numeric(1L))

        y_fit <- as.numeric(train_d[["deltaij"]])
        offset <- log(train_d[["tij"]])
        fit_meta <- list(
          basis_list = basis_list,
          encoding = encoding,
          encoded_names = colnames(X_raw)
        )

        if (.self$cross_validation) {
          prefit_args <- .self$fit_arguments
          if (!is.na(.self$maxit_prefit)) {
            prefit_args[["maxit"]] <- .self$maxit_prefit
          }

          cv_fit <- suppressWarnings(do.call(glmnet::cv.glmnet, c(
            list(
              x = X_hal,
              y = y_fit,
              offset = offset,
              penalty.factor = penalty_factor
            ),
            prefit_args
          )))

          if (length(cv_fit$lambda.min) == 0L || !is.finite(cv_fit$lambda.min)) {
            return(make_failed_fit("cv.glmnet failed in Learner_hal5::private_fit"))
          }
          .self$lambda_opt <- as.numeric(cv_fit$lambda.min)
          attach_psl_fit_meta(cv_fit, fit_meta)
        } else {
          glmnet_args <- .self$fit_arguments
          glmnet_args[["penalty.factor"]] <- penalty_factor
          fit_obj <- suppressWarnings(do.call(glmnet::glmnet, c(
            glmnet_args,
            list(x = X_hal, y = y_fit, offset = offset)
          )))
          fit_obj$selected_lambda <- tail(fit_obj$lambda, 1L)
          attach_psl_fit_meta(fit_obj, fit_meta)
        }
      }, error = function(e) {
        make_failed_fit(conditionMessage(e))
      })

      fit
    },

    private_fit_all_causes = function(data, causes, grid_nodes, ...) {
      fits <- lapply(causes, function(cause) {
        .self$private_fit(data, cause = cause, grid_nodes = grid_nodes, ...)
      })
      names(fits) <- as.character(causes)
      fits
    },

    private_predictor = function(model, newdata, grid_nodes, ...) {
      pred_cols <- unique(c(.self$basic_covariates, "node"))
      x <- data.table::copy(newdata[, .SD, .SDcols = pred_cols])
      data.table::setnames(x, "node", "terminal_node")
      x[, terminal_node_ix := match(terminal_node, grid_nodes)]

      if (anyNA(x$terminal_node_ix)) {
        stop("Some terminal_node values are not exactly in grid_nodes.")
      }

      n_expanded <- sum(x$terminal_node_ix)
      if (n_expanded == 0L) return(numeric(0L))
      if (is_failed_fit(model)) return(rep(NA_real_, n_expanded))

      fit_meta <- get_psl_fit_meta(model)
      if (is.null(fit_meta)) return(rep(NA_real_, n_expanded))

      idx <- rep.int(seq_len(nrow(x)), times = x$terminal_node_ix)
      pred_grid <- x[idx]
      pred_grid[, node_ix := sequence(x$terminal_node_ix)]
      pred_grid[, node := factor(
        node_ix,
        levels = seq_along(grid_nodes),
        labels = paste0("n", seq_along(grid_nodes))
      )]

      encoded <- list()
      for (variable in names(fit_meta$encoding)) {
        value <- pred_grid[[variable]]
        variable_encoding <- fit_meta$encoding[[variable]]
        if (identical(variable_encoding$type, "factor")) {
          value <- factor(value, levels = variable_encoding$levels)
          if (anyNA(value)) return(rep(NA_real_, n_expanded))
          if (length(variable_encoding$levels) > 1L) {
            for (level in variable_encoding$levels[-1L]) {
              column_name <- paste0(variable, "_", make.names(level))
              encoded[[column_name]] <- as.numeric(value == level)
            }
          }
        } else {
          encoded[[variable]] <- as.numeric(value)
        }
      }

      X_raw <- do.call(cbind, encoded)
      if (is.null(dim(X_raw))) X_raw <- matrix(X_raw, ncol = 1L)
      storage.mode(X_raw) <- "double"
      X_raw <- X_raw[, fit_meta$encoded_names, drop = FALSE]

      X_new <- try(
        hal9001::make_design_matrix(X_raw, fit_meta$basis_list),
        silent = TRUE
      )
      if (inherits(X_new, "try-error")) return(rep(NA_real_, n_expanded))

      pred_args <- list(
        object = model,
        newx = X_new,
        newoffset = rep(0.0, nrow(X_new)),
        type = "response"
      )
      dots <- list(...)
      if (inherits(model, "cv.glmnet") && is.null(dots[["s"]])) {
        pred_args[["s"]] <- "lambda.min"
      }
      if (!inherits(model, "cv.glmnet") && is.null(dots[["s"]]) &&
          length(model$selected_lambda) == 1L && is.finite(model$selected_lambda)) {
        pred_args[["s"]] <- model$selected_lambda
      }

      out <- try(do.call(stats::predict, c(pred_args, dots)), silent = TRUE)
      if (inherits(out, "try-error")) return(rep(NA_real_, n_expanded))
      out <- as.numeric(out)
      if (length(out) != n_expanded) return(rep(NA_real_, n_expanded))
      out
    }
  )
)
