#' \code{xgboost} learner class
#'
#' @import data.table
#' @export Learner_xgboost
#' @exportClass Learner_xgboost
Learner_xgboost <- setRefClass(
  "Learner_xgboost",
  fields = list(
    covariates = "character",
    treatment = "character",
    nrounds = "integer",
    cross_validation = "logical",
    intercept = "logical",
    formula = "character",
    learner = "function",
    add_nodes = "logical",
    fit_arguments = "list",
    cv_hyperparameters_grid = "data.frame",
    grid_of_hyperparameters = "list",
    nrounds_cv = "integer",
    nfold_cv = "integer",
    fit_arguments_cv = "list"
  ),
  methods = list(
    initialize = function(covariates = NA_character_,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          intercept = FALSE,
                          add_nodes = TRUE,
                          nrounds = NA_integer_,
                          nfold_cv = 5L,
                          nrounds_cv = 100L,
                          grid_of_hyperparameters = NULL,
                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$cross_validation <- cross_validation

      .self$add_nodes <- add_nodes

      .self$intercept <- intercept

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula(
        covariates = .self$covariates,
        treatment = .self$treatment,
        intercept = .self$intercept,
        add_nodes = .self$add_nodes
      )


      # handle fit arguments


      .self$nrounds <- nrounds

      .self$fit_arguments <- list(...)

      .self$fit_arguments <- list(params = list(objective = "count:poisson"))
      .self$fit_arguments[['nrounds']] <-  .self$nrounds

      # handle hyperparamters grid if necessary

      if (.self$cross_validation) {
        .self$fit_arguments_cv <- list(params = list(objective = "count:poisson"))
        .self$fit_arguments_cv[['nrounds']] <-  nrounds_cv
        .self$fit_arguments_cv[['nfold']] <-  nfold_cv
        .self$fit_arguments_cv[['verbose']] <-  FALSE

        .self$grid_of_hyperparameters <- grid_of_hyperparameters
        .self$cv_hyperparameters_grid  <- expand.grid(grid_of_hyperparameters)

      } else{
        .self$grid_of_hyperparameters <- grid_of_hyperparameters
        #
        .self$fit_arguments[['params']] <- c(.self$fit_arguments[['params']], grid_of_hyperparameters)

      }





    },

    update_cross_validation_argument = function(nfold) {
      .self$fit_arguments_cv[['nfold']] <- nfold

    },

    private_fit = function(data, ...) {
      x = sparse.model.matrix(formula(.self$formula), data)


      xgb.mx <- xgb.DMatrix(data = x, label = data[['deltaij']])

      setinfo(xgb.mx, "base_margin", log(data[['tij']]))

      .self$fit_arguments[['data']] <- xgb.mx



      if (.self$cross_validation) {
        .self$fit_arguments_cv[['data']] <- xgb.mx

        best_nll <- Inf

        for (i in 1:nrow(.self$cv_hyperparameters_grid)) {
          .self$fit_arguments_cv[['params']] <- c(as.list(.self$cv_hyperparameters_grid[i, ]))

          out_cv <- do.call(xgb.cv, .self$fit_arguments_cv)

          current_nll <- min(out_cv$evaluation_log$test_poisson_nloglik_mean)

          if (current_nll < best_nll) {
            best_nll <- current_nll
            best_params <- .self$fit_arguments_cv[['params']]
          }
        }


        .self$fit_arguments[['params']][names(best_params)] <- NULL
        .self$fit_arguments[['params']] <- c(.self$fit_arguments[['params']], best_params)

      }

      out <- do.call(xgb.train, .self$fit_arguments)


      return(out)

    },

    predictor = function(model, newdata, ...) {
      x = sparse.model.matrix(formula(.self$formula), newdata)

      # xgb.mx <- xgb.DMatrix(data = x,
      #                       label = newdata[['deltaij']])

      xgtest = xgb.DMatrix(x)
      setinfo(xgtest, "base_margin", log(newdata[['tij']]))

      out <- predict(model, newdata = xgtest)

      return(out)


    },
    private_predictor = function(model, newdata, ...) {
      x = sparse.model.matrix(formula(.self$formula), newdata)

      # xgb.mx <- xgb.DMatrix(data = x,
      #                       label = newdata[['deltaij']])

      xgtest = xgb.DMatrix(x)
      setinfo(xgtest, "base_margin", rep(1, nrow(newdata)))

      out <- predict(model, newdata = xgtest)

      return(out)


    }




  )
)

#' \code{glmnet} learner class
#'
#' @export Learner_glmnet
#' @exportClass Learner_glmnet
Learner_glmnet <- setRefClass(
  "Learner_glmnet",
  fields = list(
    covariates = "character",
    treatment = "character",
    cross_validation = "logical",
    recycle_information = "logical",
    intercept = "logical",
    formula = "character",
    learner = "function",
    add_nodes = "logical",
    penalise_nodes = "logical",
    lambda_grid = "numeric",
    lambda = "numeric",
    fit_arguments = "list",
    covariates_attributes_matrix = "list",
    id = "character",
    model_fit = "ANY",
    data_info = "list"
  ),
  methods = list(
    initialize = function(covariates = NA_character_,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          intercept = TRUE,
                          add_nodes = TRUE,
                          penalise_nodes = FALSE,
                          recycle_information = FALSE,
                          id = NA_character_,
                          lambda_grid = NA_real_,
                          lambda = NA_real_,
                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$id <- id

      .self$cross_validation <- cross_validation

      .self$intercept <- intercept

      .self$add_nodes <- add_nodes

      .self$penalise_nodes <- penalise_nodes

      tmp <- lambda_grid
      # normalize user input
      if (is.null(lambda_grid))
        tmp <- NA_real_
      if (length(lambda_grid) == 0L)
        tmp <- NA_real_
      .self$lambda_grid <- as.numeric(tmp)  # will keep NA as NA

      tmp <- lambda
      if (is.null(lambda))
        tmp <- NA_real_
      if (length(lambda) == 0L)
        tmp <- NA_real_
      tmp <- as.numeric(tmp)

      .self$recycle_information <- recycle_information

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula_glmnet(
        covariates = .self$covariates,
        treatment = .self$treatment,
        add_nodes = .self$add_nodes
      )



      .self$fit_arguments <- list(...)

      .self$covariates_attributes_matrix <- list(...)

      .self$fit_arguments[['family']] <- "poisson"

      .self$fit_arguments[['intercept']] <- .self$intercept

      if (.self$cross_validation) {
        .self$learner = cv.glmnet

      } else{
        .self$learner = glmnet
        .self$fit_arguments[['lambda']] <- tmp
      }


    },



    update_cross_validation_argument = function(nfold) {
      .self$fit_arguments[['nfolds']] <- nfold

    },

    save_meta_data = function(id, ...) {
      .self$id <- id

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


      group_cols <- c(.self$covariates, .self$treatment)[complete.cases(c(.self$covariates, .self$treatment))]


      group_cols <- group_cols[is.character(group_cols) &
                                 !is.na(group_cols) & nzchar(group_cols)]
      group_cols <- unlist(lapply(group_cols, .extract_symbols), use.names = FALSE)

      group_cols <- unique(group_cols)

      data <- data[, .(tij = sum(tij), deltaij = sum(deltaij)), by = c(group_cols, "node", "k")]

      data <- data[complete.cases(data), ]

      x = sparse.model.matrix(formula(.self$formula), data, contrasts.arg = NULL)[, -1]


      if (!.self$penalise_nodes) {
        .self$fit_arguments[['penalty.factor']] <- 1 - (grepl("node", colnames(x)))#& !grepl("node:", colnames(x))& !grepl(":node", colnames(x)))

      }





      if (.self$cross_validation) {
        suppressWarnings(cv_fit <- do.call(cv.glmnet, c(
          .self$fit_arguments,
          list(
            x = x,
            y = as.numeric(data[["deltaij"]]),
            offset = log(data[["tij"]])
          )
        )))

        lambda_grid_prefit <- cv_fit$lambda


        if (!all(is.na(.self$lambda_grid))) {
          lambda_grid_prefit <- c(.self$lambda_grid, lambda_grid_prefit)
          lambda_grid_prefit <- lambda_grid_prefit[complete.cases(lambda_grid_prefit)]
          lambda_grid_prefit <- sort(unique(lambda_grid_prefit), decreasing = TRUE)
        }
        preds <- glmnet::predict.glmnet(
          cv_fit$glmnet.fit,
          newx = x,
          newoffset = log(data[['tij']]),
          type = "response",
          s = lambda_grid_prefit
        )

        if (is.null(dim(preds))) {
          preds <- matrix(preds, ncol = 1L)
        }

        mu <- pmax(preds, 1e-12)
        y_vec <- as.numeric(data[['deltaij']])
        nll <- -colSums(y_vec * log(mu) - mu)
        lambda_opt <- lambda_grid_prefit[which.min(nll)]

        glmnet_args <- .self$fit_arguments
        glmnet_args[['lambda']] <- lambda_opt
        glmnet_args[['nfolds']] <- NULL

        out <- do.call(glmnet, c(glmnet_args, list(
          x = x,
          y = as.numeric(data[["deltaij"]]),
          offset = log(data[["tij"]])
        )))

      } else {
        out <- do.call(.self$learner, c(
          .self$fit_arguments,
          list(
            x = x,
            y = as.numeric(data[["deltaij"]]),
            offset = log(data[["tij"]])
          )
        ))
      }

      return(out)

    },

    predictor = function(model, newdata, ...) {
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
        ...,
        newx = sparse.model.matrix(formula(.self$formula), newdata, contrasts.arg = NULL)[, -1],
        newoffset = log(newdata[['tij']]),
        type = "response"
      )

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
        newx = sparse.model.matrix(formula(.self$formula), newdata, contrasts.arg = NULL)[, -1],
        newoffset = log(1),
        type = "response",
        ...
      )

      return(out)


    }


  #   predict = function(newdata, times, cause = 1, ...) {
  #
  #     setDT(newdata)
  #
  #     tmp <- data.table::copy(newdata)
  #     # here we disregard the event_time column if present in the newdata
  #     tmp <- tmp[, setdiff(
  #       names(tmp),
  #       c(
  #         .self$data_info$event_time,
  #         .self$data_info$status,
  #         .self$data_info$id
  #       )
  #     ), with = FALSE]
  #
  #
  #     ## checks on the data
  #     cond_zero <- 0 %in% times
  #
  #     all_zero <- all(times == 0)
  #
  #     cond_times_larger_than_max <- times > .self$data_info$maximum_followup
  #
  #     ## frame hazard problem
  #     pwch_cols <- paste0("pwch_", 1:.self$data_info$n_crisks)
  #
  #     if (all(cond_times_larger_than_max)) {
  #       warning(
  #         paste0(
  #           "All the entries in the input times are larger than the maximum follow-up: ",
  #           as.character(.self$data_info$maximum_followup)
  #         )
  #       )
  #       d <- NULL
  #
  #       return(d)
  #
  #     } else{
  #       eval(parse(
  #         text = paste0(
  #           "
  #   vec_dt <- data.table(
  #
  #   ",
  #           .self$data_info$event_time,
  #           " = times[times <= .self$data_info$maximum_followup]
  # )
  #   "
  #         )
  #       ))
  #
  #       # # no problem writing over id
  #       # if (is.null(tmp[[.self$data_info$id]])) {
  #       #   tmp[[.self$data_info$id]] <- 1:nrow(tmp)
  #       # }
  #       #
  #       # if (is.null(tmp[[.self$data_info$status]])) {
  #       #   tmp[[.self$data_info$status]] <- 0
  #       # }
  #
  #       tmp[, dummy := 1]
  #       vec_dt[, dummy := 1]
  #
  #       # Merge on dummy to create Cartesian product
  #       data_pp <- merge(tmp, vec_dt, by = "dummy", allow.cartesian = TRUE)[, dummy := NULL]
  #     }
  #     # }
  #
  #
  #
  #     if (cond_zero) {
  #       zero_time <- data_pp[time == 0, ]
  #       zero_time[, c(pwch_cols, 'survival_function', 'absolute_risk') := list(rep(0, length(pwch_cols)), 1, 0)]
  #       data_pp <- data_pp[time != 0, ]
  #     }
  #
  #     if (all_zero) {
  #       return(zero_time)
  #     }
  #
  #     # no problem writing over id
  #     if (is.null(data_pp[[.self$data_info$id]])) {
  #       data_pp[[.self$data_info$id]] <- 1:nrow(data_pp)
  #     }
  #
  #     if (is.null(data_pp[[.self$data_info$status]])) {
  #       data_pp[[.self$data_info$status]] <- 0
  #     }
  #
  #     data_pp <- data_pre_processing(
  #       data_pp,
  #       id = .self$data_info$id,
  #       status = .self$data_info$status,
  #       predictions = TRUE,
  #       event_time = .self$data_info$event_time,
  #       nodes = .self$data_info$nodes
  #     )
  #
  #
  #
  #     #  ()
  #     if (!is.null(.self$data_info$variable_transformation)) {
  #
  #       # Take the variable that we transform
  #       apply_transformations(data_pp, .self$data_info$variable_transformation)
  #
  #
  #     }
  #
  #
  #     # Predict on the validation set your pseudo-observations ----
  #     #
  #
  #     data_pp[, deltatime := tij][, tij := 1]
  #
  #
  #     dt_pred <-lapply(.self$model_fit,
  #                      function(x) {
  #                        .self$private_predictor(model=x,
  #                                                newdata = data_pp)
  #                      })
  #
  #     browser()
  #     # dt_pred <- lapply(object$superlearner, function(sl_fit) {
  #     #   object$learners[[1]]$private_predictor(model = sl_fit$learners_fit, newdata = data_pp)
  #     # })
  #
  #
  #     # save casue-specific pwch
  #
  #     data_pp[, paste0("pwch_", 1:.self$data_info$n_crisks) := dt_pred]
  #
  #
  #
  #     # save sum of pwch
  #
  #     sum_of_hazards <- paste(pwch_cols, collapse = " + ")
  #
  #     pwch_dot_string <- paste0("data_pp[, pwch_dot :=", sum_of_hazards, "]")
  #
  #     eval(parse(text = pwch_dot_string))
  #
  #
  #     # compute cumulative hazard
  #
  #     mapply(function(pwch, name) {
  #       data_pp[, (paste0("cumulative_hazard_", name)) := cumsum(get(pwch) * deltatime), by = id]
  #     }, pwch_cols, gsub("pwch_", "", pwch_cols))
  #
  #
  #     # compute survival function
  #
  #     ## c++
  #     haz <- as.matrix(data_pp[, .SD, .SDcols = patterns("^pwch_[0-9]+$")])
  #     S <- pch_survival(id = data_pp$id,
  #                       dt = data_pp$deltatime,
  #                       haz = haz)
  #     data_pp[, survival_function := S]
  #
  #
  #     data_pp[, absolute_risk := pch_absolute_risk(id, deltatime, haz, cause_idx = cause)]
  #
  #     # abs_risk approx
  #     # data_pp[, survival_function_shift := shift(survival_function, fill = 1), by =
  #     #           id]
  #     # absolute_risk_string <- paste0(
  #     #   "data_pp[, absolute_risk_2 := cumsum(survival_function_shift * pwch_",
  #     #   cause,
  #     #   "*deltatime), by = id]"
  #     # )
  #     # eval(parse(text = absolute_risk_string))
  #
  #
  #     ## non c++
  #     # hazard_terms <- paste0("cumulative_hazard_", 1:.self$data_info$n_crisks)
  #     # sum_expr <- paste(pwch_cols, collapse = " + ")
  #     # survival_function_string <- paste0("data_pp[, survival_function := exp(-cumsum((", sum_expr, ")*deltatime)),by=id]")
  #     # eval(parse(text = survival_function_string))
  #
  #     # shift survival function
  #     # data_pp[, survival_function_shift := shift(survival_function, fill = 1), by =
  #     #           id]
  #     # absolute_risk_string <- paste0(
  #     #   "data_pp[, absolute_risk := cumsum(survival_function_shift * pwch_",
  #     #   cause,
  #     #   "/pwch_dot * (1-exp(-pwch_dot*deltatime))), by = id]"
  #     # )
  #     # eval(parse(text = absolute_risk_string))
  #
  #     ####
  #     data_pp <- data_pp[, .SD[.N], by = id][, times := as.numeric(as.character(node)) +
  #                                              deltatime]
  #
  #
  #
  #     columns_ss <- unique(
  #       c(
  #         colnames(newdata),
  #         .self$data_info$event_time,
  #         pwch_cols,
  #         "survival_function",
  #         "absolute_risk"
  #       )
  #     )
  #
  #     d <- data_pp[, ..columns_ss]
  #
  #
  #     if (cond_zero) {
  #       d <- rbind(zero_time, d)
  #
  #     }
  #
  #
  #
  #     if (any(cond_times_larger_than_max)) {
  #       eval(parse(
  #         text = paste0(
  #           "
  #   vec_dt2 <- data.table(
  #
  #   ",
  #           .self$data_info$event_time,
  #           " = times[cond_times_larger_than_max]
  # )
  #   "
  #         )
  #       ))
  #
  #       tmp[, dummy := 1]
  #       vec_dt2[, dummy := 1]
  #
  #       d2 <- merge(tmp, vec_dt2, by = "dummy", allow.cartesian = TRUE)[, dummy := NULL]
  #       d2[, c(pwch_cols, 'survival_function') := list(rep(NA, length(pwc_cols)), NA)]
  #
  #
  #       if (.self$data_info$id %in% colnames(d)) {
  #         d2[[.self$data_info$id]] <- (nrow(data_pp) + 1):(nrow(data_pp) + nrow(d2))
  #       }
  #
  #
  #
  #       d <- rbind(d, d2)
  #
  #
  #     }
  #
  #
  #     return(d)
  #
  #
  #
  #
  #   }
  #




  )
)

#' \code{hal} learner class
#'
#' @export Learner_hal
#' @exportClass Learner_hal
Learner_hal <- setRefClass(
  "Learner_hal",
  fields = list(
    covariates = "character",
    treatment = "character",
    cross_validation = "logical",
    recycle_information = "logical",
    intercept = "logical",
    formula = "character",
    learner = "function",
    add_nodes = "logical",
    max_degree = "integer",
    lambda_grid = "numeric",
    lambda_opt = "numeric",
    penalise_nodes = "logical",
    maxit_prefit = "numeric",
    fit_arguments = "list",
    covariates_attributes_matrix = "list",
    num_knots = "numeric",
    id = "character",
    model_fit = "ANY"
  ),
  methods = list(
    initialize = function(covariates = NA_character_,
                          treatment = NA_character_,
                          id = NA_character_,
                          intercept = FALSE,
                          cross_validation = TRUE,
                          add_nodes = TRUE,
                          penalise_nodes = FALSE,
                          recycle_information = FALSE,
                          max_degree = 2,
                          maxit_prefit = NA_real_,
                          num_knots = c(10L, 5L),
                          lambda_grid = NA_real_,
                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$id <- id

      .self$cross_validation <- cross_validation

      .self$intercept <- intercept

      .self$add_nodes <- add_nodes

      .self$max_degree <- max_degree

      .self$num_knots <- num_knots

      # normalize user input
      tmp <- lambda_grid
      if (is.null(lambda_grid))
        tmp <- NA_real_
      if (length(lambda_grid) == 0L)
        tmp <- NA_real_

      .self$lambda_grid <- tmp

      .self$penalise_nodes <- penalise_nodes

      .self$recycle_information <- recycle_information

      .self$maxit_prefit <- maxit_prefit

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula_hal(
        covariates = .self$covariates,
        treatment = .self$treatment,
        intercept = FALSE,
        #in the glmnet case, intercept is handled separately.
        add_nodes = .self$add_nodes
      )

      if (.self$cross_validation) {
        .self$learner <- cv.glmnet
      } else {
        .self$learner <- glmnet
      }

      .self$fit_arguments <- list(...)

      .self$covariates_attributes_matrix <- list()

      .self$fit_arguments[['family']] <- "poisson"

      .self$fit_arguments[['intercept']] <- .self$intercept




    },

    save_meta_data = function(id, ...) {
      .self$id <- id

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


    update_cross_validation_argument = function(nfold) {
      .self$fit_arguments[['nfolds']] <- nfold

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

      group_cols <- c(.self$covariates, .self$treatment)[complete.cases(c(.self$covariates, .self$treatment))]

      data_copy <- data_copy[, .(tij = sum(tij), deltaij = sum(deltaij)), by = c(group_cols, "node", "k")]

      data_copy <- data_copy[complete.cases(data_copy), ]

      x_pp <- hal_basis(
        vars = c(group_cols, "node"),
        DT = data_copy,
        max_interaction = .self$max_degree,
        knots_per_order = .self$num_knots
      )

      if (!.self$penalise_nodes) {
        .self$fit_arguments[['penalty.factor']] <- 1 - grepl('^I\\(\\s*node\\s*==[^)]*\\)$', x_pp$colnames)

      }

      .self$covariates_attributes_matrix <- x_pp[c("colnames", "meta")]

      #
      # if (.self$cross_validation) {
      #
      #   prefit_args <- .self$fit_arguments
      #
      #   if (!is.na(.self$maxit_prefit)) {
      #     prefit_args[['maxit']] <- .self$maxit_prefit
      #   }
      #
      #   suppressWarnings(cv_fit <- do.call(cv.glmnet, c(prefit_args,
      #                                  list(y=as.numeric(data_copy[['deltaij']]),
      #                                       offset=log(data_copy[['tij']]),
      #                                       x=x_pp$X,
      #                                       maxit=5000
      #                                       ))))
      #
      #   lambda_grid_prefit <- cv_fit$lambda
      #
      #   if (!all(is.na(.self$lambda_grid))) {
      #     lambda_grid_prefit <- c(.self$lambda_grid, lambda_grid_prefit)
      #     lambda_grid_prefit <- lambda_grid_prefit[complete.cases(lambda_grid_prefit)]
      #     lambda_grid_prefit <- sort(unique(lambda_grid_prefit), decreasing = TRUE)
      #   }
      #
      #   preds <- predict(cv_fit$glmnet.fit,
      #                    newx = x_pp$X,
      #                    newoffset = log(data_copy[['tij']]),
      #                    type = "response",
      #                    s = lambda_grid_prefit)
      #
      #   if (is.null(dim(preds))) {
      #     preds <- matrix(preds, ncol = 1L)
      #   }
      #
      #   mu <- pmax(preds, 1e-12)
      #   y_vec <- as.numeric(data_copy[['deltaij']])
      #   nll <- -colSums(y_vec * log(mu) - mu)
      #   # lambda_opt <- lambda_grid_prefit[which.min(nll)]
      #
      #   glmnet_args <- .self$fit_arguments
      #   glmnet_args[['lambda']] <- lambda_opt
      #   glmnet_args[['nfolds']] <- NULL
      #
      #   fit <- cv_fit
      #   .self$lambda_opt <- lambda_grid_prefit[which.min(nll)]
      #   fit <- do.call(glmnet, c(glmnet_args,
      #                            list(y=as.numeric(data_copy[['deltaij']]),
      #                                 offset=log(data_copy[['tij']]),
      #                                 x=x_pp$X
      #                            )))
      # }
      #

      if (.self$cross_validation) {
        ## 1) Prefit with cv.glmnet to obtain the lambda path + CV curve
        prefit_args <- .self$fit_arguments

        if (!is.na(.self$maxit_prefit)) {
          prefit_args[["maxit"]] <- .self$maxit_prefit
        }

        suppressWarnings(cv_fit <- do.call(cv.glmnet, c(
          prefit_args, list(
            x      = x_pp$X,
            y      = as.numeric(data_copy[["deltaij"]]),
            offset = log(data_copy[["tij"]])
          )
        )))

        ## 2) Choose lambda that minimizes cvm (this is cv_fit$lambda.min)
        lambda_min_cvm <- cv_fit$lambda.min
        .self$lambda_opt <- lambda_min_cvm

        ## 3) Refit a plain glmnet at that single lambda (no CV)
        glmnet_args <- .self$fit_arguments
        glmnet_args[["lambda"]] <- lambda_min_cvm
        glmnet_args[["nfolds"]] <- NULL  # not used by glmnet, but remove if present

        suppressWarnings(fit <- do.call(glmnet, c(
          glmnet_args, list(
            x      = x_pp$X,
            y      = as.numeric(data_copy[["deltaij"]]),
            offset = log(data_copy[["tij"]])
          )
        )))

      } else {
        fit <- tryCatch(
          withCallingHandlers(
            do.call(.self$learner, .self$fit_arguments),
            warning = function(w) {
              if (grepl("Convergence", conditionMessage(w), ignore.case = TRUE) ||
                  grepl("empty model", conditionMessage(w), ignore.case = TRUE)) {
                invokeRestart("muffleWarning")
              }
            }
          ),
          error = function(e) {
            return(NULL)
          }
        )
      }

      return(fit)

    },


    predictor = function(model, newdata, ...) {
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
        ...,
        newx = X_new,
        newoffset = log(newdata[['tij']]),
        type = "response"
      )
      #}

      return(out)


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


#' \code{gam} learner class using mgcv bam
#'
#' @param covariates \code{character}
#' @param treatment \code{character}
#' @param cross_validation \code{logical}
#'
#'
#'
#' @export Learner_gam
#' @exportClass Learner_gam
Learner_gam <- setRefClass(
  "Learner_gam",
  fields = list(
    covariates = "character",
    treatment = "character",
    cross_validation = "logical",
    recycle_information = "logical",
    intercept = "logical",
    formula = "character",
    learner = "function",
    add_nodes = "logical",
    penalise_nodes = "logical",
    fit_arguments = "list",
    covariates_attributes_matrix = "list",
    id = "character",
    model_fit = "ANY"
  ),
  methods = list(
    initialize = function(covariates = NULL,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          intercept = TRUE,
                          add_nodes = TRUE,
                          penalise_nodes = FALSE,
                          recycle_information = FALSE,
                          id = NA_character_,
                          ...) {
      .self$covariates <- covariates
      .self$treatment <- treatment
      .self$id <- id
      .self$cross_validation <- cross_validation
      .self$intercept <- intercept
      .self$add_nodes <- add_nodes
      .self$penalise_nodes <- penalise_nodes
      .self$recycle_information <- recycle_information

      .self$formula <- create_formula_gam(
        covariates = .self$covariates,
        treatment = .self$treatment,
        intercept = .self$intercept,
        add_nodes = .self$add_nodes
      )

      .self$learner <- mgcv::bam

      .self$fit_arguments <- list(...)
      .self$fit_arguments[['family']] <- poisson()
    },

    update_cross_validation_argument = function(nfold) {

    },

    save_meta_data = function(id, ...) {
      .self$id <- id

    },

    datapp = function(train_data = NULL,
                      validation_data = NULL) {
      if (!is.null(train_data)) {
        train.mf <- model.frame(
          as.formula(.self$formula),
          rbind(train_data, validation_data),
          drop.unused.levels = FALSE
        )
        x <- model.matrix(attr(train.mf, "terms"), data = train_data)
        y <- train_data[['deltaij']]
        offset <- log(train_data[['tij']])
        .self$covariates_attributes_matrix[['train.mf']] <- train.mf
        .self$recycle_information <- TRUE
      } else {
        x <- model.matrix(attr(.self$covariates_attributes_matrix[['train.mf']], "terms"),
                          data = validation_data)
        y <- validation_data[['deltaij']]
        offset <- log(validation_data[['tij']])
      }

      return(list(x = x, y = y, offset = offset))
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


      group_cols <- c(.self$covariates, .self$treatment)[complete.cases(c(.self$covariates, .self$treatment))]


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

    fit = function(data,
                   id = "id",
                   #
                   stratified_k_fold = FALSE,
                   #
                   start_time = NULL,
                   #
                   end_time = NULL,
                   #
                   status = "status",
                   #
                   event_time = NULL,
                   #
                   number_of_nodes = NULL,
                   #
                   nodes = NULL,
                   variable_transformation = NULL,
                   nfold = 3,
                   ...) {
      # Multiple checks about the input ----
      ############
      check_1 <- is.null(start_time) & !is.null(end_time)
      check_2 <- !is.null(start_time) & is.null(end_time)
      check_3 <- (!is.null(start_time) ||
                    !is.null(end_time)) & !is.null(event_time)


      if (!(id %in% names(data))) {
        data[["id"]] <- 1:NROW(data)
        id <<- "id"
      }



      if (check_1 || check_2) {
        stop("For interval data, both start_time and end_time are required.")

      }

      if (check_3) {
        stop("Either provide interval data or censored data")

      }

      if (!is.null(start_time) & !is.null(end_time)) {
        interval_data_type = TRUE

        maximum_followup = max(data[[end_time]])

      } else{
        interval_data_type = FALSE

        maximum_followup = max(data[[event_time]])

      }

      n <- length(unique(data[[id]]))

      if (!(0 %in% data[[status]])) {
        warning(
          paste0(
            "There is no value of ",
            status,
            " equal to zero. We will consider the data uncensored."
          )
        )
        n_crisks <- length(unique(data[[status]]))
        uncensored_01 <- TRUE

      } else{
        n_crisks <- length(unique(data[[status]])) - 1
        uncensored_01 <- FALSE
      }


      if (interval_data_type) {
        # Compute here the nodes checking
        if (!is.null(number_of_nodes)) {
          observed_times <- c(data[[start_time]], data[[end_time]])
          grid_nodes <- seq(min(observed_times),
                            max(observed_times) + 1,
                            length.out = as.integer(number_of_nodes))

        } else{
          if (is.null(nodes)) {
            observed_times <- c(data[[start_time]], data[[end_time]])
            grid_nodes <- sort(unique(observed_times))
          } else {
            grid_nodes <- sort(nodes)
          }
        }

        if (!(0 %in% grid_nodes)) {
          grid_nodes <- c(0, grid_nodes)
        }

        # Actual data pp
        dt <- data_pre_processing_interval_data(
          data = data,
          id = id,
          status = status,
          start_time = start_time,
          nodes = grid_nodes,
          end_time = end_time
        )




      } else{
        #  Handle nodes
        ##Either the nodes are given or we take all of the realised times
        if (!is.null(number_of_nodes)) {
          # grid_nodes <- seq(min(data[[event_time]]), max(data[[event_time]]) + 1, length.out = as.integer(number_of_nodes))
          # grid_nodes <- unique(sort(sample(data[[event_time]],
          #                      as.integer(number_of_nodes))))

          grid_nodes = quantile(
            data[[event_time]],
            probs = seq(0, 1, length.out = as.integer(number_of_nodes) + 1),
            type = 1,
            names = FALSE
          )



        } else{
          if (is.null(nodes)) {
            grid_nodes <- sort(unique(data[[event_time]]))

          } else{
            grid_nodes <- nodes

          }
        }

        # Add zero if missing
        if (!(0 %in% grid_nodes)) {
          grid_nodes <- c(0, grid_nodes)

        }



        grid_nodes <- grid_nodes[grid_nodes <= max(data[[event_time]])]

        # Actual data pp
        dt <- data_pre_processing(
          data = data,
          id = id,
          status = status,
          nodes = grid_nodes,
          event_time = event_time,
          uncensored_01 = uncensored_01
        )



      }


      lhs_string = NULL

      if (!is.null(variable_transformation)) {
        apply_transformations(dt, variable_transformation)
      }


      # .self$model_fit <- .self$private_fit(dt)

      training_data <- split(dt, by = "k")

      .self$model_fit <- mapply(function(x) {
        out <- .self$private_fit(x)
        return(out)
      }, training_data, SIMPLIFY = FALSE)


    },

    predictor = function(model, newdata, ...) {
      pred <- predict(
        model,
        newdata = newdata,
        type = "response",
        offset = log(newdata[['tij']]),
        ...
      )
      return(pred)
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
