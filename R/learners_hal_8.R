#' HAL learner using the PICASSO Poisson solver
#'
#' `Learner_hal8` retains the piecewise-constant hazard preprocessing, HAL basis
#' construction, and prediction interface of [Learner_hal()], but fits the
#' resulting Poisson model with [picasso::picasso()] or
#' [picasso::cv.picasso()]. PICASSO requires a dense numeric design matrix, so
#' the sparse HAL basis is converted immediately before fitting and prediction.
#'
#' @details
#' PICASSO does not provide `glmnet`'s `penalty.factor` argument. Consequently,
#' all HAL basis columns, including node main effects, receive the same L1
#' penalty in `Learner_hal8`.
#'
#' @inheritParams Learner_hal
#' @return An initialized `Learner_hal8` Reference Class object.
#' @export Learner_hal8
#' @exportClass Learner_hal8
Learner_hal8 <- methods::setRefClass(
  "Learner_hal8",
  contains = "Learner_hal",
  methods = list(
    initialize = function(...) {
      callSuper(...)
      if (!requireNamespace("picasso", quietly = TRUE)) {
        stop("Learner_hal8 requires the 'picasso' package.")
      }
      .self$learner <- if (.self$cross_validation) {
        picasso::cv.picasso
      } else {
        picasso::picasso
      }
    },

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

        x_pp <- .self$hal_basis(
          vars = c(.self$basic_covariates, "node"),
          DT = train_d,
          max_interaction = .self$max_degree,
          knots_per_order = .self$num_knots
        )
        if (nrow(x_pp$X) == 0L || ncol(x_pp$X) == 0L) {
          return(make_failed_fit("Empty HAL basis in Learner_hal8::private_fit"))
        }

        X_hal <- as.matrix(x_pp$X)
        storage.mode(X_hal) <- "double"
        y <- as.numeric(train_d[["deltaij"]])
        offset <- log(train_d[["tij"]])
        fit_meta <- x_pp[c("colnames", "meta")]

        fit_args <- hal8_picasso_arguments(
          .self$fit_arguments,
          maxit_prefit = .self$maxit_prefit,
          cross_validation = .self$cross_validation
        )
        base_args <- list(
          X = X_hal,
          Y = y,
          family = "poisson",
          method = "l1",
          intercept = .self$intercept,
          offset = offset
        )

        fit_obj <- suppressWarnings(do.call(
          if (.self$cross_validation) picasso::cv.picasso else picasso::picasso,
          c(base_args, fit_args)
        ))

        if (.self$cross_validation) {
          if (length(fit_obj$lambda.min) != 1L ||
              !is.finite(fit_obj$lambda.min)) {
            return(make_failed_fit(
              "cv.picasso failed in Learner_hal8::private_fit"
            ))
          }
          .self$lambda_opt <- as.numeric(fit_obj$lambda.min)
        } else {
          fit_obj$selected_lambda <- tail(fit_obj$lambda, 1L)
        }

        attach_psl_fit_meta(fit_obj, fit_meta)
      }, error = function(e) {
        make_failed_fit(conditionMessage(e))
      })
      fit
    },

    private_fit_all_causes = function(data, causes, grid_nodes, ...) {
      tryCatch({
        prepared <- hal8_prepare_all_causes(
          learner = .self,
          data = data,
          causes = causes,
          grid_nodes = grid_nodes
        )
        X_hal <- as.matrix(prepared$X)
        storage.mode(X_hal) <- "double"
        fit_args <- hal8_picasso_arguments(
          .self$fit_arguments,
          maxit_prefit = .self$maxit_prefit,
          cross_validation = .self$cross_validation
        )
        base_args <- list(
          X = X_hal,
          family = "poisson",
          method = "l1",
          intercept = .self$intercept,
          offset = prepared$offset
        )

        fits <- lapply(seq_along(causes), function(index) {
          fit_obj <- suppressWarnings(do.call(
            if (.self$cross_validation) {
              picasso::cv.picasso
            } else {
              picasso::picasso
            },
            c(base_args, list(Y = prepared$y[, index]), fit_args)
          ))

          if (.self$cross_validation) {
            if (length(fit_obj$lambda.min) != 1L ||
                !is.finite(fit_obj$lambda.min)) {
              return(make_failed_fit(paste0(
                "cv.picasso failed in Learner_hal8::private_fit_all_causes for cause ",
                causes[index]
              )))
            }
            .self$lambda_opt <- as.numeric(fit_obj$lambda.min)
          } else {
            fit_obj$selected_lambda <- tail(fit_obj$lambda, 1L)
          }
          attach_psl_fit_meta(fit_obj, prepared$fit_meta)
        })
        names(fits) <- as.character(causes)
        fits
      }, error = function(e) {
        stats::setNames(
          rep(list(make_failed_fit(conditionMessage(e))), length(causes)),
          as.character(causes)
        )
      })
    },

    private_predictor = function(model, newdata, grid_nodes, ...) {
      pred_cols <- unique(c(.self$basic_covariates, "node"))
      x <- data.table::copy(newdata[, .SD, .SDcols = pred_cols])
      data.table::setnames(x, "node", "terminal_node")
      x[, terminal_node_ix := match(terminal_node, grid_nodes)]

      if (anyNA(x$terminal_node_ix)) {
        stop("Some terminal_node values are not exactly in grid_nodes.")
      }
      if (any(x$terminal_node_ix < 1L)) {
        stop("Some terminal_node_ix values are smaller than 1.")
      }

      n_expanded <- sum(x$terminal_node_ix)
      if (n_expanded == 0L) return(numeric(0L))
      if (is_failed_fit(model)) return(rep(NA_real_, n_expanded))

      fit_meta <- get_psl_fit_meta(model)
      if (is.null(fit_meta)) return(rep(NA_real_, n_expanded))

      idx <- rep.int(seq_len(nrow(x)), times = x$terminal_node_ix)
      pred_grid <- x[idx]
      pred_grid[, node_ix := sequence(x$terminal_node_ix)]
      pred_grid[, `:=`(tij = 1.0, deltaij = 0.0)]
      pred_grid[, node := factor(
        node_ix,
        levels = seq_along(grid_nodes),
        labels = paste0("n", seq_along(grid_nodes))
      )]

      X_new <- try(.self$hal_prepare_new(pred_grid, fit_meta), silent = TRUE)
      if (inherits(X_new, "try-error") || is.null(X_new) ||
          nrow(X_new) == 0L || ncol(X_new) == 0L) {
        return(rep(NA_real_, n_expanded))
      }
      X_new <- as.matrix(X_new)
      storage.mode(X_new) <- "double"

      dots <- list(...)
      pred_args <- list(
        object = model,
        newdata = X_new,
        type = "response",
        newoffset = rep(0, nrow(X_new))
      )
      if (inherits(model, "cv.picasso") && is.null(dots[["s"]])) {
        pred_args[["s"]] <- "lambda.min"
      }
      if (!inherits(model, "cv.picasso") && is.null(dots[["s"]]) &&
          length(model$selected_lambda) == 1L &&
          is.finite(model$selected_lambda)) {
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

hal8_picasso_arguments <- function(arguments, maxit_prefit,
                                   cross_validation) {
  arguments[["family"]] <- NULL
  arguments[["intercept"]] <- NULL
  arguments[["method"]] <- NULL

  if (!is.null(arguments[["alpha"]])) {
    if (length(arguments$alpha) != 1L || !isTRUE(arguments$alpha == 1)) {
      stop("Learner_hal8 supports Lasso (`alpha = 1`) only.")
    }
    arguments[["alpha"]] <- NULL
  }
  if (!is.null(arguments[["penalty.factor"]])) {
    stop("PICASSO does not support `penalty.factor`.")
  }
  if (!is.null(arguments[["thresh"]]) && is.null(arguments[["prec"]])) {
    arguments[["prec"]] <- arguments[["thresh"]]
  }
  arguments[["thresh"]] <- NULL
  if (!is.null(arguments[["maxit"]]) && is.null(arguments[["max.ite"]])) {
    arguments[["max.ite"]] <- arguments[["maxit"]]
  }
  arguments[["maxit"]] <- NULL

  if (isTRUE(cross_validation) && !is.na(maxit_prefit) &&
      is.null(arguments[["max.ite"]])) {
    arguments[["max.ite"]] <- maxit_prefit
  }
  arguments
}

# This is Learner_hal's shared all-cause preprocessing, factored out only so the
# final solver call can be switched from glmnet to PICASSO.
hal8_prepare_all_causes <- function(learner, data, causes, grid_nodes) {
  widths <- c(diff(grid_nodes), 1.0)
  needed_cols <- unique(c(learner$basic_covariates, "node", "tij", "deltaij"))
  train_d <- data.table::copy(data[, .SD, .SDcols = needed_cols])

  if (length(learner$basic_covariates) > 0L) {
    group_key <- unique(train_d[, .SD, .SDcols = learner$basic_covariates])
    group_key[, gid := .I]
    train_d <- group_key[train_d, on = learner$basic_covariates]
  } else {
    group_key <- data.table::data.table(gid = 1L)
    train_d[, gid := 1L]
  }
  train_d <- train_d[!is.na(gid) & !is.na(node) & !is.na(tij)]
  train_d[, node_ix := match(node, grid_nodes)]

  event_counts <- train_d[
    deltaij %in% causes,
    .N,
    by = .(gid, node_ix, deltaij)
  ]
  terminal_base <- train_d[, .(
    tij = sum(tij),
    N_terminal = .N
  ), by = .(gid, node, node_ix)]
  data.table::setorder(terminal_base, gid, node_ix)
  terminal_base[, deltaij := 0.0]

  expanded <- expand_terminal_grouped_cpp(
    gid = terminal_base$gid,
    node_ix = terminal_base$node_ix,
    tij = terminal_base$tij,
    deltaij = terminal_base$deltaij,
    N = terminal_base$N_terminal,
    widths = widths
  )
  train_long <- data.table::as.data.table(expanded)
  train_long[, node := grid_nodes[node_ix]]
  train_long <- group_key[train_long, on = "gid"]
  data.table::setorder(train_long, gid, node_ix)
  train_long[, node := factor(
    node_ix,
    levels = seq_along(grid_nodes),
    labels = paste0("n", seq_along(grid_nodes))
  )]

  x_pp <- learner$hal_basis(
    vars = c(learner$basic_covariates, "node"),
    DT = train_long,
    max_interaction = learner$max_degree,
    knots_per_order = learner$num_knots
  )
  if (nrow(x_pp$X) == 0L || ncol(x_pp$X) == 0L) {
    stop("Empty HAL basis in Learner_hal8::private_fit_all_causes")
  }

  cause_names <- as.character(causes)
  train_index <- train_long[, .(gid, node_ix)]
  if (nrow(event_counts) > 0L) {
    event_wide <- data.table::dcast(
      event_counts,
      gid + node_ix ~ deltaij,
      value.var = "N",
      fill = 0L
    )
  } else {
    event_wide <- unique(train_index)
  }
  missing_causes <- setdiff(cause_names, names(event_wide))
  if (length(missing_causes) > 0L) {
    event_wide[, (missing_causes) := 0L]
  }
  event_wide <- event_wide[train_index, on = .(gid, node_ix)]
  for (cause_name in cause_names) {
    missing <- which(is.na(event_wide[[cause_name]]))
    if (length(missing) > 0L) {
      data.table::set(event_wide, missing, cause_name, 0L)
    }
  }

  list(
    X = x_pp$X,
    y = as.matrix(event_wide[, ..cause_names]),
    offset = log(train_long[["tij"]]),
    fit_meta = x_pp[c("colnames", "meta")]
  )
}
