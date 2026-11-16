#' HAL learner fitted directly with `hal9001::fit_hal()`
#'
#' `Learner_hal6` uses the same piecewise-constant Poisson hazard data as
#' [Learner_hal()], but passes those data directly to [hal9001::fit_hal()].
#' Thus, `hal9001` performs both zero-order basis construction and model fitting;
#' no precomputed HAL design matrix is supplied by this package.
#'
#' @inheritParams Learner_hal
#' @return An initialized `Learner_hal6` Reference Class object.
#' @export Learner_hal6
#' @exportClass Learner_hal6
Learner_hal6 <- methods::setRefClass(
  "Learner_hal6",
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

        expanded <- expand_terminal_grouped_cpp(
          gid = train_d$gid,
          node_ix = train_d$node_ix,
          tij = train_d$tij,
          deltaij = train_d$deltaij,
          N = train_d$N_terminal,
          widths = widths
        )

        train_d <- data.table::as.data.table(expanded)
        train_d <- group_key[train_d, on = "gid"]
        data.table::setorder(train_d, gid, node_ix)
        train_d[, node := factor(
          node_ix,
          levels = seq_along(grid_nodes),
          labels = paste0("n", seq_along(grid_nodes))
        )]

        encoded <- hal6_encode_training_data(
          train_d,
          variables = c(.self$basic_covariates, "node")
        )
        x <- encoded$x
        y <- as.numeric(train_d[["deltaij"]])

        if (nrow(x) == 0L || ncol(x) == 0L) {
          return(make_failed_fit("Empty input data in Learner_hal6::private_fit"))
        }

        fit_args <- .self$fit_arguments
        lambda <- fit_args[["lambda"]]
        fit_args[c("family", "lambda", "offset", "weights")] <- NULL
        fit_args[["cv_select"]] <- .self$cross_validation
        # HAL's default response bounds are based on event counts. They are not
        # valid bounds for hazards predicted with an offset of zero.
        fit_args[["prediction_bounds"]] <- c(-Inf, Inf)
        if (!is.na(.self$maxit_prefit) && .self$cross_validation) {
          fit_args[["maxit"]] <- .self$maxit_prefit
        }

        hal_fit_0 <- suppressWarnings(hal9001::fit_hal(
          X = x,
          Y = y,
          max_degree = .self$max_degree,
          smoothness_orders = 0,
          num_knots = .self$num_knots,
          family = "poisson",
          lambda = lambda,
          offset = log(train_d[["tij"]]),
          fit_control = fit_args
        ))

        # Without CV, fit_hal retains every point on the glmnet path. The base
        # learner interface requires one prediction per row, so use the least
        # penalized (last) fitted model, as Learner_hal does.
        if (!.self$cross_validation && ncol(hal_fit_0$coefs) > 1L) {
          selected <- ncol(hal_fit_0$coefs)
          hal_fit_0$coefs <- hal_fit_0$coefs[, selected, drop = FALSE]
          hal_fit_0$lambda_star <- tail(hal_fit_0$lambda_star, 1L)
        }

        attach_psl_fit_meta(
          hal_fit_0,
          list(
            encoding = encoded$encoding,
            encoded_names = colnames(x)
          )
        )
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

      new_x <- try(
        hal6_encode_new_data(
          pred_grid,
          encoding = fit_meta$encoding,
          encoded_names = fit_meta$encoded_names
        ),
        silent = TRUE
      )
      if (inherits(new_x, "try-error")) return(rep(NA_real_, n_expanded))

      out <- try(
        stats::predict(
          model,
          new_data = new_x,
          offset = rep(0, nrow(new_x)),
          type = "response",
          ...
        ),
        silent = TRUE
      )
      if (inherits(out, "try-error")) return(rep(NA_real_, n_expanded))

      out <- as.numeric(out)
      if (length(out) != n_expanded) return(rep(NA_real_, n_expanded))
      out
    }
  )
)

# Convert factor-like variables to a stable numeric input matrix for fit_hal.
hal6_encode_training_data <- function(data, variables) {
  encoding <- setNames(vector("list", length(variables)), variables)
  encoded <- list()

  for (variable in variables) {
    value <- data[[variable]]
    if (is.character(value) || is.logical(value)) value <- factor(value)

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

  x <- do.call(cbind, encoded)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  storage.mode(x) <- "double"
  list(x = x, encoding = encoding)
}

hal6_encode_new_data <- function(data, encoding, encoded_names) {
  encoded <- list()

  for (variable in names(encoding)) {
    value <- data[[variable]]
    variable_encoding <- encoding[[variable]]
    if (identical(variable_encoding$type, "factor")) {
      value <- factor(value, levels = variable_encoding$levels)
      if (anyNA(value)) stop("New data contain an unseen factor level.")
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

  x <- do.call(cbind, encoded)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  storage.mode(x) <- "double"
  x[, encoded_names, drop = FALSE]
}
