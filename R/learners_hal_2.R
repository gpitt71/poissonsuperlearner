#' HAL learner using custom sparse Poisson IRLS
#'
#' `Learner_hal_2` reuses [Learner_hal()] basis construction and prediction
#' metadata, but fits its Lasso path with [hal_poisson_irls()] instead of
#' `glmnet`. This Phase 1 implementation operates on the existing `dgCMatrix`
#' HAL design and does not use token or cumulative-threshold representations.
#'
#' @inheritParams Learner_hal
#' @param ... Solver arguments including `lambda`, `tol`, `maxit_irls`,
#'   `maxit_cd`, `active_set`, `nlambda`, and `lambda.min.ratio`.
#' @return An initialized `Learner_hal_2` Reference Class object.
#' @export Learner_hal_2
#' @exportClass Learner_hal_2
Learner_hal_2 <- methods::setRefClass(
  "Learner_hal_2",
  contains = "Learner_hal",
  methods = list(
    initialize = function(covariates = NA_character_,
                          intercept = TRUE,
                          cross_validation = TRUE,
                          max_degree = 2L,
                          maxit_prefit = NA_real_,
                          num_knots = c(10L, 5L),
                          ...) {
      callSuper(covariates = covariates, intercept = intercept,
                cross_validation = cross_validation, max_degree = max_degree,
                maxit_prefit = maxit_prefit, num_knots = num_knots, ...)
      .self$learner <- if (cross_validation) cv_hal_poisson_irls else hal_poisson_irls
      .self$fit_arguments$family <- NULL
      if (!is.null(.self$fit_arguments$alpha) && .self$fit_arguments$alpha != 1)
        stop("Learner_hal_2 supports Lasso (`alpha = 1`) only.")
      if (isTRUE(.self$fit_arguments$standardize))
        stop("Learner_hal_2 Phase 1 requires `standardize = FALSE`.")
      .self$fit_arguments$alpha <- NULL
      .self$fit_arguments$standardize <- NULL
    },

    private_fit = function(data, cause, grid_nodes, ...) {
      .self$private_fit_all_causes(data, causes = cause,
                                   grid_nodes = grid_nodes, ...)[[1L]]
    },

    private_fit_all_causes = function(data, causes, grid_nodes, ...) {
      tryCatch({
        prep <- hal2_prepare_training_data(
          data = data, basic_covariates = .self$basic_covariates,
          causes = causes, grid_nodes = grid_nodes)
        basis <- .self$hal_basis(
          vars = c(.self$basic_covariates, "node"), DT = prep$data,
          max_interaction = .self$max_degree,
          knots_per_order = .self$num_knots)
        if (!nrow(basis$X) || !ncol(basis$X))
          stop("Empty HAL basis in Learner_hal_2.")

        penalty_factor <- 1 - grepl(
          "^I\\(\\s*node\\s*==[^)]*\\)$", basis$colnames)
        args <- .self$fit_arguments
        args$intercept <- .self$intercept
        args$penalty_factor <- penalty_factor
        if (!is.null(args$lambda.min.ratio)) {
          args$lambda_min_ratio <- args$lambda.min.ratio
          args$lambda.min.ratio <- NULL
        }
        if (!is.na(.self$maxit_prefit) && is.null(args$maxit_cd))
          args$maxit_cd <- as.integer(.self$maxit_prefit)
        if (.self$cross_validation && is.null(args$nfolds)) args$nfolds <- 5L

        fit_meta <- basis[c("colnames", "meta")]
        output <- setNames(vector("list", length(causes)), as.character(causes))
        for (index in seq_along(causes)) {
          output[[index]] <- tryCatch({
            fit <- do.call(.self$learner, c(list(
              X = basis$X,
              y = as.numeric(prep$y[, index]),
              offset = log(prep$data$tij)), args))
            if (.self$cross_validation) .self$lambda_opt <- fit$lambda.min
            if (!.self$cross_validation)
              fit$selected_lambda <- tail(fit$lambda, 1L)
            attach_psl_fit_meta(fit, fit_meta)
          }, error = function(e) {
            make_failed_fit(paste0("Learner_hal_2 failed for cause ", causes[index],
                                   ": ", conditionMessage(e)))
          })
        }
        output
      }, error = function(e) {
        setNames(rep(list(make_failed_fit(conditionMessage(e))), length(causes)),
                 as.character(causes))
      })
    }
  )
)

hal2_prepare_training_data <- function(data, basic_covariates, causes, grid_nodes) {
  needed <- unique(c(basic_covariates, "node", "tij", "deltaij"))
  terminal <- data.table::copy(data[, .SD, .SDcols = needed])
  if (length(basic_covariates)) {
    key <- unique(terminal[, .SD, .SDcols = basic_covariates])
    key[, gid := .I]
    terminal <- key[terminal, on = basic_covariates]
  } else {
    key <- data.table::data.table(gid = 1L)
    terminal[, gid := 1L]
  }
  terminal <- terminal[!is.na(gid) & !is.na(node) & is.finite(tij)]
  terminal[, node_ix := match(node, grid_nodes)]
  if (anyNA(terminal$node_ix)) stop("Training node is not in `grid_nodes`.")
  events <- terminal[deltaij %in% causes, .N, by = .(gid, node_ix, deltaij)]
  skeleton <- terminal[, .(tij = sum(tij), N_terminal = .N),
                       by = .(gid, node, node_ix)]
  data.table::setorder(skeleton, gid, node_ix)
  skeleton[, deltaij := 0]
  expanded <- data.table::as.data.table(expand_terminal_grouped_cpp(
    gid = skeleton$gid, node_ix = skeleton$node_ix, tij = skeleton$tij,
    deltaij = skeleton$deltaij, N = skeleton$N_terminal,
    widths = c(diff(grid_nodes), 1)))
  expanded <- key[expanded, on = "gid"]
  data.table::setorder(expanded, gid, node_ix)
  expanded[, node := factor(node_ix, levels = seq_along(grid_nodes),
                            labels = paste0("n", seq_along(grid_nodes)))]
  if (any(!is.finite(expanded$tij) | expanded$tij <= 0))
    stop("Learner_hal_2 requires finite, positive aggregated exposure.")

  index <- expanded[, .(gid, node_ix)]
  wide <- if (nrow(events)) {
    data.table::dcast(events, gid + node_ix ~ deltaij, value.var = "N", fill = 0L)
  } else unique(index)
  cause_names <- as.character(causes)
  missing <- setdiff(cause_names, names(wide))
  if (length(missing)) wide[, (missing) := 0L]
  wide <- wide[index, on = .(gid, node_ix)]
  y <- as.matrix(wide[, ..cause_names])
  y[is.na(y)] <- 0
  list(data = expanded, y = y)
}
