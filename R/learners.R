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
    intercept ="logical",
    formula ="character",
    learner="function",
    add_nodes="logical",
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
                          add_nodes= TRUE,
                          nrounds = NA_integer_,
                          nfold_cv = 5L,
                          nrounds_cv = 100L,
                          grid_of_hyperparameters=NULL,
                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$cross_validation <- cross_validation

      .self$add_nodes <- add_nodes

      .self$intercept<-intercept

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula(covariates = .self$covariates,
                                      treatment = .self$treatment,
                                      intercept = .self$intercept,
                                      add_nodes=.self$add_nodes)


      # handle fit arguments


      .self$nrounds<-nrounds

      .self$fit_arguments <- list(...)

      .self$fit_arguments <- list(params=list(objective="count:poisson"))
      .self$fit_arguments[['nrounds']] <-  .self$nrounds

      # handle hyperparamters grid if necessary

      if (.self$cross_validation) {

        .self$fit_arguments_cv <- list(params=list(objective="count:poisson"))
        .self$fit_arguments_cv[['nrounds']] <-  nrounds_cv
        .self$fit_arguments_cv[['nfold']] <-  nfold_cv
        .self$fit_arguments_cv[['verbose']] <-  FALSE

        .self$grid_of_hyperparameters <- grid_of_hyperparameters
        .self$cv_hyperparameters_grid  <- expand.grid(grid_of_hyperparameters)

      } else{

        .self$grid_of_hyperparameters <-grid_of_hyperparameters
        #
        .self$fit_arguments[['params']] <- c(.self$fit_arguments[['params']],
                                             grid_of_hyperparameters)

      }





    },

    update_cross_validation_argument= function(nfold){

      .self$fit_arguments_cv[['nfold']] <- nfold

    },

    fit = function(data, ...) {

      x = sparse.model.matrix(formula(.self$formula),
                              data)


      xgb.mx <- xgb.DMatrix(data = x,
                            label = data[['deltaij']])

      setinfo(xgb.mx,
              "base_margin",
              log(data[['tij']]))

      .self$fit_arguments[['data']] <-xgb.mx



      if( .self$cross_validation){

        .self$fit_arguments_cv[['data']] <-xgb.mx

        best_nll <- Inf

        for(i in 1:nrow(.self$cv_hyperparameters_grid)){

          .self$fit_arguments_cv[['params']] <- c(as.list(.self$cv_hyperparameters_grid[i,]))

          out_cv <- do.call(xgb.cv,
                     .self$fit_arguments_cv)

          current_nll <- min(out_cv$evaluation_log$test_poisson_nloglik_mean)

          if (current_nll < best_nll) {
            best_nll <- current_nll
            best_params <- .self$fit_arguments_cv[['params']]
          }
        }


        .self$fit_arguments[['params']][names(best_params)]<-NULL
        .self$fit_arguments[['params']] <-c(.self$fit_arguments[['params']],
                                            best_params)

      }

      out <- do.call(xgb.train,
                     .self$fit_arguments)


      return(out)

    },

    predictor = function(model, newdata, ...) {


      x = sparse.model.matrix(formula(.self$formula),
                              newdata)

      # xgb.mx <- xgb.DMatrix(data = x,
      #                       label = newdata[['deltaij']])

      xgtest = xgb.DMatrix(x)
      setinfo(xgtest, "base_margin", log(newdata[['tij']]))

      out <- predict(model,
                     newdata = xgtest)

      return(out)


    },
    private_predictor = function(model, newdata, ...) {


      x = sparse.model.matrix(formula(.self$formula),
                              newdata)

      # xgb.mx <- xgb.DMatrix(data = x,
      #                       label = newdata[['deltaij']])

      xgtest = xgb.DMatrix(x)
      setinfo(xgtest, "base_margin", rep(1,nrow(newdata)))

      out <- predict(model,
                     newdata = xgtest)

      return(out)


    }




    ))

#' \code{glm} learner class
#'
#' @export Learner_glm
#' @exportClass Learner_glm
Learner_glm <- setRefClass(
  "Learner_glm",
  fields = list(
    covariates = "character",
    treatment = "character",
    cross_validation = "logical",
    intercept ="logical",
    formula ="character",
    learner="function",
    add_nodes="logical"
  ),
  methods = list(
    initialize = function(covariates = NA_character_,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          intercept = TRUE,
                          add_nodes= TRUE,
                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$cross_validation <- cross_validation

      .self$add_nodes <- add_nodes

      .self$intercept<-intercept

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula(covariates = .self$covariates,
                                      treatment = .self$treatment,
                                      intercept = .self$intercept,
                                      add_nodes=.self$add_nodes)


      if (.self$cross_validation) {
        .self$learner = glm
        warning("There is no cross-validation procedure for a glm model. This learner is a simple glm.")

      } else{
        .self$learner = glm
      }
    },

    fit = function(data, ...) {


      out <- glm.fit(x=sparse.model.matrix(formula(.self$formula),
                                           data),#[ss,]
                     y=data[['deltaij']], #[ss]
                     offset = log(data[['tij']]), #[ss]
                     family = poisson())

      return(out)

    },

    # fit = function(data, validation_data=NULL, ...) {

      # survival_01 <- length(unique(data[['k']])) < 2


      # practical correction in case it is a survival problem
      # if(survival_01){
      #
      #   .self$formula <- create_formula(covariates = .self$covariates,
      #                                   treatment = .self$treatment,
      #                                   competing_risks =FALSE,
      #                                   add_nodes=.self$add_nodes)
      #
      # }




    #   out<- .self$learner(.self$formula,
    #                       data=rbind(data,
    #                                  validation_data),
    #                       family = "poisson",
    #                       subset = rep(TRUE,dim(data)[1]))
    #
    #   return(out)
    #
    # },


    fit = function(data, ...) {

      # out<- .self$learner(.self$formula,
      #                     data=data,
      #                     family = "poisson", ...)

      # ss <- data[['train_01']]
      out <- glm.fit(x=sparse.model.matrix(formula(.self$formula),
                                           data),#[ss,]
                     y=data[['deltaij']], #[ss]
                     offset = log(data[['tij']]), #[ss]
                     family = poisson())

      return(out)

    },

    predictor = function(model, newdata,...) {

      newdata<-newdata[complete.cases(newdata),]

      newx=sparse.model.matrix(formula(.self$formula),
                               newdata)

      newoffset = log(newdata[['tij']])

      out <-exp(newx %*% model$coefficients + newoffset)


      return(as.numeric(out))



    },

    private_predictor = function(model, newdata, ...) {



      newdata<-newdata[complete.cases(newdata),]


      newx=sparse.model.matrix(formula(.self$formula),
                               newdata)#[ss,]

      newoffset = log(newdata[['tij']])#[ss])


      out <-exp(newx %*% model$coefficients + newoffset)


      return(as.numeric(out))






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
    intercept="logical",
    formula ="character",
    learner="function",
    add_nodes="logical",
    penalise_nodes= "logical",
    lambda_grid="numeric",
    fit_arguments = "list",
    covariates_attributes_matrix= "list",
    id="character"
  ),
  methods = list(

    initialize = function(covariates = NA_character_,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          intercept=FALSE,
                          add_nodes = TRUE,
                          penalise_nodes=FALSE,
                          recycle_information =FALSE,
                          id=NA_character_,
                          lambda_grid=NA_real_,
                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$id <- id

      .self$cross_validation <- cross_validation

      .self$intercept <- intercept

      .self$add_nodes <- add_nodes

      .self$penalise_nodes <- penalise_nodes

      .self$lambda_grid <- lambda_grid

      .self$recycle_information <- recycle_information

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula(covariates = .self$covariates,
                                      treatment = .self$treatment,
                                      intercept =FALSE, #in the glmnet case, intercept is handled separately.
                                      add_nodes=.self$add_nodes)


      if (.self$cross_validation) {
        .self$learner = cv.glmnet

      } else{
        .self$learner = glmnet
      }

      .self$fit_arguments <- list(...)

      .self$covariates_attributes_matrix <- list(...)

      .self$fit_arguments[['family']] <- "poisson"

      .self$fit_arguments[['intercept']] <- .self$intercept




    },



    update_cross_validation_argument= function(nfold){

      .self$fit_arguments[['nfolds']] <- nfold

    },

    save_meta_data= function(id, ...){

      .self$id <- id

    },

    fit = function(data, ...) {

      .extract_symbols <- function(term) {
        # try as a bare expression, then as a RHS of a formula
        out <- tryCatch(
          all.vars(str2lang(term)),
          error = function(e) {
            tryCatch(
              all.vars(stats::terms(stats::as.formula(paste("~", term)))),
              error = function(e2) character(0)
            )
          }
        )
        out
      }


      group_cols <- c(.self$covariates,.self$treatment)[complete.cases(c(.self$covariates,.self$treatment))]


      group_cols <- group_cols[is.character(group_cols) & !is.na(group_cols) & nzchar(group_cols)]
      group_cols <- unlist(lapply(group_cols, .extract_symbols), use.names = FALSE)

      group_cols <- unique(group_cols)

      data <- data[, .(tij = sum(tij), deltaij = sum(deltaij)), by = c(group_cols, "node", "k")]

      data<-data[complete.cases(data),]

      x = sparse.model.matrix(formula(.self$formula),
                              data,
                              contrasts.arg = NULL)

      .self$fit_arguments[['y']] <- as.numeric(data[['deltaij']])

      .self$fit_arguments[['offset']] <-  log(data[['tij']])


      if(!.self$penalise_nodes){

        .self$fit_arguments[['penalty.factor']] <- 1- (grepl("node",colnames(x)))#& !grepl("node:", colnames(x))& !grepl(":node", colnames(x)))

      }

      .self$fit_arguments[['x']] <- x

      if (.self$cross_validation) {
        cv_fit <- do.call(cv.glmnet, .self$fit_arguments)

        if (is.null(cv_fit$glmnet.fit)) return(cv_fit)

        lambda_grid <- cv_fit$lambda

        if (!all(is.na(.self$lambda_grid))) {
          lambda_grid <- c(.self$lambda_grid, lambda_grid)
          lambda_grid <- lambda_grid[complete.cases(lambda_grid)]
          lambda_grid <- sort(unique(lambda_grid), decreasing = TRUE)
        }
        preds <- predict(cv_fit$glmnet.fit,
                         newx = x,
                         newoffset = log(data[['tij']]),
                         type = "response",
                         s = lambda_grid)

        if (is.null(dim(preds))) {
          preds <- matrix(preds, ncol = 1L)
        }

        mu <- pmax(preds, 1e-12)
        y_vec <- as.numeric(data[['deltaij']])
        nll <- -colSums(y_vec * log(mu) - mu)
        lambda_opt <- lambda_grid[which.min(nll)]

        glmnet_args <- .self$fit_arguments
        glmnet_args[['lambda']] <- lambda_opt
        glmnet_args[['nfolds']] <- NULL

        out <- do.call(glmnet, glmnet_args)

      } else {
        out <- do.call(.self$learner,
                       .self$fit_arguments)
      }

      return(out)

    },


    predictor = function(model, newdata, ...) {

      is_empty_model <- function(model) {
        if (is.null(model)) return(TRUE)
        if (inherits(model, "cv.glmnet")) {
          if (is.null(model$glmnet.fit)) return(TRUE)
          beta <- model$glmnet.fit$beta
        } else {
          beta <- model$beta
        }
        if (is.null(beta)) return(TRUE)
        if (!is.matrix(beta) && !inherits(beta, "Matrix")) return(FALSE)
        nrow(beta) == 0 || ncol(beta) == 0
      }

      if (is_empty_model(model)) {
        return(rep(NA_real_, nrow(newdata)))
      }

      out <- predict(model,
                     ...,
                   newx=sparse.model.matrix(formula(.self$formula),
                                            newdata),
                   newoffset = log(newdata[['tij']]),
                   type = "response")

      return(out)


    },

    private_predictor = function(model, newdata, ...) {


      is_empty_model <- function(model) {
        if (is.null(model)) return(TRUE)
        if (inherits(model, "cv.glmnet")) {
          if (is.null(model$glmnet.fit)) return(TRUE)
          beta <- model$glmnet.fit$beta
        } else {
          beta <- model$beta
        }
        if (is.null(beta)) return(TRUE)
        if (!is.matrix(beta) && !inherits(beta, "Matrix")) return(FALSE)
        nrow(beta) == 0 || ncol(beta) == 0
      }

      if (is_empty_model(model)) {
        return(rep(NA_real_, nrow(newdata)))
      }

      out <- predict(model,
                     ...,
                     newx=sparse.model.matrix(formula(.self$formula),
                                              newdata),#[ss,],
                     newoffset = log(1),#log(newdata[['tij']]),#[ss]),
                     type = "response")


      # if(all((levels(newdata$node) %in% model$xlevels$node))){
      #
      #   return(pred)
      #
      # }else{
      #
      #
      #   #
      #   newdata[node %in% model$xlevels$node,predictions_model:=pred]
      #
      #   return(as.array(newdata$predictions_model))
      #
      #
      #
      #
      # }

      return(out)


    }



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
    intercept="logical",
    formula ="character",
    learner="function",
    add_nodes="logical",
    max_degree="integer",
    lambda_grid="numeric",
    penalise_nodes= "logical",
    maxit_prefit="numeric",
    fit_arguments = "list",
    covariates_attributes_matrix= "list",
    num_knots="numeric",
    id="character"
  ),
  methods = list(

    initialize = function(covariates = NA_character_,
                          treatment = NA_character_,
                          id = NA_character_,
                          intercept=FALSE,
                          cross_validation=TRUE,
                          add_nodes = TRUE,
                          penalise_nodes=FALSE,
                          recycle_information =FALSE,
                          max_degree=2,
                          maxit_prefit=NA_real_,
                          num_knots=c(10L,5L),
                          lambda_grid=NA_real_,
                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$id <- id

      .self$cross_validation <- cross_validation

      .self$intercept <- intercept

      .self$add_nodes <- add_nodes

      .self$max_degree <- max_degree

      .self$num_knots <- num_knots

      .self$lambda_grid <- lambda_grid

      .self$penalise_nodes <- penalise_nodes

      .self$recycle_information <- recycle_information

      .self$maxit_prefit <- maxit_prefit

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula_hal(covariates = .self$covariates,
                                          treatment = .self$treatment,
                                          intercept =FALSE, #in the glmnet case, intercept is handled separately.
                                          add_nodes=.self$add_nodes)

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

    save_meta_data = function(id, ...){

      .self$id <- id

    },

    hal_basis = function(vars,
                          DT,
                          max_interaction = 2L,
                          knots_per_order = rep(10L, 2L)) {

      stopifnot(requireNamespace("data.table", quietly = TRUE))
      stopifnot(requireNamespace("Matrix", quietly = TRUE))

      if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)
      if (length(knots_per_order) != max_interaction)
        stop("knots_per_order must have length max_interaction")
      if (!all(vars %in% names(DT))) stop("vars not all found in DT")

      n <- nrow(DT)
      DT_local <- DT[, ..vars]
      for (v in vars) {
        x <- DT_local[[v]]
        if (is.character(x) || is.logical(x)) DT_local[[v]] <- factor(x)
      }

      inter2 <- function(a, b) {
        i <- 1L; j <- 1L; la <- length(a); lb <- length(b)
        if (!la || !lb) return(integer())
        out <- integer(min(la, lb)); k <- 0L
        while (i <= la && j <= lb) {
          ai <- a[i]; bj <- b[j]
          if (ai < bj) i <- i + 1L else if (ai > bj) j <- j + 1L else {
            k <- k + 1L; out[k] <- ai; i <- i + 1L; j <- j + 1L
          }
        }
        if (k) out[seq_len(k)] else integer()
      }
      interN <- function(lst) {
        if (length(lst) == 1L) return(lst[[1L]])
        acc <- inter2(lst[[1L]], lst[[2L]])
        if (!length(acc) || length(lst) == 2L) return(acc)
        for (i in 3L:length(lst)) {
          acc <- inter2(acc, lst[[i]])
          if (!length(acc)) break
        }
        acc
      }

      primitives <- setNames(vector("list", length(vars)), vars)
      meta_per_var <- setNames(vector("list", length(vars)), vars)

      mk_main <- function(x, v, K) {
        if (is.numeric(x)) {
          ps  <- seq(0, 1, length.out = K + 2L)
          cps <- unique(stats::quantile(x, probs = ps, na.rm = TRUE, type = 1))
          cps <- cps[is.finite(cps)]
          if (length(cps) >= 2L) cps <- cps[-c(1L, length(cps))] else cps <- numeric()
          cps <- unique(cps)
          if (!length(cps)) {
            return(list(idxs = list(), names = character(),
                        prim_meta = list(), var_meta = list(cutpoints = numeric())))
          }
          ord <- order(x); xs <- x[ord]
          pos <- findInterval(cps, xs, rightmost.closed = TRUE, all.inside = TRUE)
          idxs <- lapply(pos, function(k) if (k > 0L) ord[seq_len(k)] else integer())
          idxs <- lapply(idxs, sort.int, method = "quick")
          nms  <- sprintf("I(%s<=%.10g)", v, cps)
          prim_meta <- lapply(cps, function(cu) list(var = v, kind = "numeric", cutpoint = cu))
          list(idxs = idxs, names = nms,
               prim_meta = prim_meta,
               var_meta = list(cutpoints = cps))
        } else if (is.factor(x)) {
          lv <- levels(x)
          if (!length(lv)) {
            return(list(idxs = list(), names = character(),
                        prim_meta = list(), var_meta = list(levels = character())))
          }
          # TAKE ALL LEVELS, ignore K
          idx_by_lvl <- split(seq_len(length(x)), x, drop = FALSE)
          idxs <- lapply(lv, function(l) sort.int(idx_by_lvl[[l]], method = "quick"))
          nms  <- sprintf("I(%s==%s)", v, make.names(lv, unique = TRUE))
          prim_meta <- lapply(lv, function(L) list(var = v, kind = "factor", level = L))
          list(idxs = idxs, names = nms,
               prim_meta = prim_meta,
               var_meta = list(levels = lv))
        } else {
          stop(sprintf("Unsupported type for %s", v))
        }
      }

      for (v in vars) {
        res <- mk_main(DT_local[[v]], v, K = knots_per_order[1L])
        primitives[[v]] <- res
        meta_per_var[[v]] <- res$var_meta
      }

      I_chunks <- vector("list", 1024L)
      J_chunks <- vector("list", 1024L)
      name_chunks <- vector("list", 1024L)
      colmeta_chunks <- vector("list", 1024L)
      chunk_used <- 0L
      p <- 0L

      add_cols <- function(idxs_list, nm_vec, col_meta_list) {
        keep <- lengths(idxs_list) > 0L
        if (!any(keep)) return()
        idxs_list <- idxs_list[keep]; nm_vec <- nm_vec[keep]; col_meta_list <- col_meta_list[keep]
        kcols <- length(idxs_list)
        need <- chunk_used + kcols
        if (need > length(I_chunks)) {
          newlen <- max(need, length(I_chunks) * 2L)
          length(I_chunks) <<- newlen
          length(J_chunks) <<- newlen
          length(name_chunks) <<- newlen
          length(colmeta_chunks) <<- newlen
        }
        for (i in seq_len(kcols)) {
          p <<- p + 1L
          chunk_used <<- chunk_used + 1L
          I_chunks[[chunk_used]] <<- idxs_list[[i]]
          J_chunks[[chunk_used]] <<- rep.int(p, length(idxs_list[[i]]))
          name_chunks[[chunk_used]] <<- nm_vec[[i]]
          colmeta_chunks[[chunk_used]] <<- col_meta_list[[i]]
        }
      }

      # mains
      for (v in vars) {
        prim <- primitives[[v]]
        if (!length(prim$idxs)) next
        mains_meta <- lapply(prim$prim_meta, function(m) c(list(type = "main"), m))
        add_cols(prim$idxs, prim$names, mains_meta)
      }

      # interactions: cap numerics by Kcap, include all factor levels
      if (max_interaction >= 2L) {
        for (d in 2:max_interaction) {
          if (length(vars) < d) break
          cmbs <- utils::combn(vars, d, simplify = FALSE)
          Kcap <- knots_per_order[d]
          for (cmb in cmbs) {
            prim_d <- lapply(cmb, function(v) {
              idxs <- primitives[[v]]$idxs
              pm   <- primitives[[v]]$prim_meta
              if (!length(idxs)) return(list(idxs = list(), names = character(), meta = list()))
              # decide cap per variable
              is_num <- length(pm) > 0L && identical(pm[[1]]$kind, "numeric")
              ncap <- if (is_num) min(length(idxs), Kcap) else length(idxs)  # factors: all levels
              list(idxs = idxs[seq_len(ncap)],
                   names = primitives[[v]]$names[seq_len(ncap)],
                   meta  = pm[seq_len(ncap)])
            })
            if (any(vapply(prim_d, function(z) length(z$idxs), integer(1L)) == 0L)) next

            lens <- vapply(prim_d, function(z) length(z$idxs), integer(1L))
            total <- prod(lens); if (!total) next

            for (t in 0:(total - 1L)) {
              sel <- integer(d); rem <- t
              for (j in seq_len(d)) { sel[j] <- (rem %% lens[j]) + 1L; rem <- rem %/% lens[j] }
              chosen_idx <- vector("list", d)
              nm_parts <- character(d)
              parts_meta <- vector("list", d)
              for (j in seq_len(d)) {
                chosen_idx[[j]] <- prim_d[[j]]$idxs[[sel[j]]]
                nm_parts[j] <- prim_d[[j]]$names[[sel[j]]]
                parts_meta[[j]] <- prim_d[[j]]$meta[[sel[j]]]
              }
              idx_prod <- interN(chosen_idx)
              if (!length(idx_prod)) next
              add_cols(
                list(idx_prod),
                paste0(nm_parts, collapse = ":"),
                list(list(type = "interaction", order = d, parts = parts_meta))
              )
            }
          }
        }
      }

      if (chunk_used == 0L) stop("No basis columns created.")

      I_chunks <- I_chunks[seq_len(chunk_used)]
      J_chunks <- J_chunks[seq_len(chunk_used)]
      name_chunks <- name_chunks[seq_len(chunk_used)]
      colmeta_chunks <- colmeta_chunks[seq_len(chunk_used)]

      i_vec <- unlist(I_chunks, use.names = FALSE)
      j_vec <- unlist(J_chunks, use.names = FALSE)
      names_all <- unlist(name_chunks, use.names = FALSE)

      X <- Matrix::sparseMatrix(
        i = if (length(i_vec)) i_vec else integer(),
        j = if (length(j_vec)) j_vec else integer(),
        x = if (length(i_vec)) 1L else integer(),
        dims = c(n, if (length(j_vec)) max(j_vec) else 0L),
        dimnames = list(NULL, names_all)
      )

      list(
        X = X,
        colnames = names_all,
        meta = list(
          per_var = meta_per_var,
          per_col = colmeta_chunks
        )
      )
    },


    # hal_basis = function(vars,
    #                       DT,
    #                       max_interaction = 2L,
    #                       knots_per_order = rep(10L, 2L)) {
    #
    #   stopifnot(requireNamespace("data.table", quietly = TRUE))
    #   stopifnot(requireNamespace("Matrix", quietly = TRUE))
    #
    #   if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)
    #   if (length(knots_per_order) != max_interaction)
    #     stop("knots_per_order must have length max_interaction")
    #   if (!all(vars %in% names(DT))) stop("vars not all found in DT")
    #
    #   n <- nrow(DT)
    #   DT_local <- DT[, ..vars]
    #   for (v in vars) {
    #     x <- DT_local[[v]]
    #     if (is.character(x) || is.logical(x)) DT_local[[v]] <- factor(x)
    #   }
    #
    #   inter2 <- function(a, b) {
    #     i <- 1L; j <- 1L; la <- length(a); lb <- length(b)
    #     if (la == 0L || lb == 0L) return(integer())
    #     out <- integer(min(la, lb)); k <- 0L
    #     while (i <= la && j <= lb) {
    #       ai <- a[i]; bj <- b[j]
    #       if (ai < bj) { i <- i + 1L
    #       } else if (ai > bj) { j <- j + 1L
    #       } else { k <- k + 1L; out[k] <- ai; i <- i + 1L; j <- j + 1L }
    #     }
    #     if (k) out[seq_len(k)] else integer()
    #   }
    #   interN <- function(lst) {
    #     if (length(lst) == 1L) return(lst[[1L]])
    #     acc <- inter2(lst[[1L]], lst[[2L]])
    #     if (length(acc) == 0L || length(lst) == 2L) return(acc)
    #     for (i in 3L:length(lst)) {
    #       acc <- inter2(acc, lst[[i]])
    #       if (!length(acc)) break
    #     }
    #     acc
    #   }
    #
    #   primitives <- setNames(vector("list", length(vars)), vars)
    #   meta_per_var <- setNames(vector("list", length(vars)), vars)
    #
    #   mk_main <- function(x, v, K) {
    #     if (is.numeric(x)) {
    #       ps  <- seq(0, 1, length.out = K + 2L)
    #       cps <- unique(stats::quantile(x, probs = ps, na.rm = TRUE, type = 7))
    #       cps <- cps[is.finite(cps)]
    #       if (length(cps) >= 2L) cps <- cps[-c(1L, length(cps))] else cps <- numeric()
    #       cps <- unique(cps)
    #       if (!length(cps)) {
    #         return(list(idxs = list(), names = character(),
    #                     prim_meta = list(), var_meta = list(cutpoints = numeric())))
    #       }
    #       ord <- order(x); xs <- x[ord]
    #       pos <- findInterval(cps, xs, rightmost.closed = TRUE, all.inside = TRUE)
    #       idxs <- lapply(pos, function(k) if (k > 0L) ord[seq_len(k)] else integer())
    #       idxs <- lapply(idxs, sort.int, method = "quick")
    #       nms  <- sprintf("I(%s<=%.10g)", v, cps)
    #       prim_meta <- lapply(cps, function(cu) list(var = v, kind = "numeric", cutpoint = cu))
    #       list(idxs = idxs, names = nms,
    #            prim_meta = prim_meta,
    #            var_meta = list(cutpoints = cps))
    #     } else if (is.factor(x)) {
    #       lv <- levels(x)
    #       if (!length(lv)) {
    #         return(list(idxs = list(), names = character(),
    #                     prim_meta = list(), var_meta = list(levels = character())))
    #       }
    #       # keep ALL levels; cap by K if requested
    #       lv_use <- if (is.finite(K) && K < length(lv)) head(lv, K) else lv
    #       idx_by_lvl <- split(seq_len(length(x)), x, drop = FALSE)
    #       idxs <- lapply(lv_use, function(l) sort.int(idx_by_lvl[[l]], method = "quick"))
    #       nms  <- sprintf("I(%s==%s)", v, make.names(lv_use, unique = TRUE))
    #       prim_meta <- lapply(lv_use, function(lv) list(var = v, kind = "factor", level = lv))
    #       list(idxs = idxs, names = nms,
    #            prim_meta = prim_meta,
    #            var_meta = list(levels = lv_use))
    #     } else {
    #       stop(sprintf("Unsupported type for %s", v))
    #     }
    #   }
    #
    #   for (v in vars) {
    #     res <- mk_main(DT_local[[v]], v, K = knots_per_order[1L])
    #     primitives[[v]] <- res
    #     meta_per_var[[v]] <- res$var_meta
    #   }
    #
    #   I_chunks <- vector("list", 1024L)
    #   J_chunks <- vector("list", 1024L)
    #   name_chunks <- vector("list", 1024L)
    #   colmeta_chunks <- vector("list", 1024L)
    #   chunk_used <- 0L
    #   p <- 0L
    #
    #   add_cols <- function(idxs_list, nm_vec, col_meta_list) {
    #     keep <- lengths(idxs_list) > 0L
    #     if (!any(keep)) return()
    #     idxs_list <- idxs_list[keep]; nm_vec <- nm_vec[keep]; col_meta_list <- col_meta_list[keep]
    #     kcols <- length(idxs_list)
    #     need <- chunk_used + kcols
    #     if (need > length(I_chunks)) {
    #       newlen <- max(need, length(I_chunks) * 2L)
    #       length(I_chunks) <<- newlen
    #       length(J_chunks) <<- newlen
    #       length(name_chunks) <<- newlen
    #       length(colmeta_chunks) <<- newlen
    #     }
    #     for (i in seq_len(kcols)) {
    #       p <<- p + 1L
    #       chunk_used <<- chunk_used + 1L
    #       I_chunks[[chunk_used]] <<- idxs_list[[i]]
    #       J_chunks[[chunk_used]] <<- rep.int(p, length(idxs_list[[i]]))
    #       name_chunks[[chunk_used]] <<- nm_vec[[i]]
    #       colmeta_chunks[[chunk_used]] <<- col_meta_list[[i]]
    #     }
    #   }
    #
    #   for (v in vars) {
    #     prim <- primitives[[v]]
    #     if (!length(prim$idxs)) next
    #     mains_meta <- lapply(prim$prim_meta, function(m) c(list(type = "main"), m))
    #     add_cols(prim$idxs, prim$names, mains_meta)
    #   }
    #
    #   if (max_interaction >= 2L) {
    #     for (d in 2:max_interaction) {
    #       if (length(vars) < d) break
    #       cmbs <- utils::combn(vars, d, simplify = FALSE)
    #       Kcap <- knots_per_order[d]
    #       for (cmb in cmbs) {
    #         prim_d <- lapply(cmb, function(v) {
    #           idxs <- primitives[[v]]$idxs
    #           pm   <- primitives[[v]]$prim_meta
    #           if (!length(idxs)) return(list(idxs = list(), names = character(), meta = list()))
    #           ncap <- min(length(idxs), Kcap)
    #           list(idxs = idxs[seq_len(ncap)],
    #                names = primitives[[v]]$names[seq_len(ncap)],
    #                meta  = pm[seq_len(ncap)])
    #         })
    #         if (any(vapply(prim_d, function(z) length(z$idxs), integer(1L)) == 0L)) next
    #
    #         lens <- vapply(prim_d, function(z) length(z$idxs), integer(1L))
    #         total <- prod(lens); if (!total) next
    #
    #         for (t in 0:(total - 1L)) {
    #           sel <- integer(d); rem <- t
    #           for (j in seq_len(d)) { sel[j] <- (rem %% lens[j]) + 1L; rem <- rem %/% lens[j] }
    #           chosen_idx <- vector("list", d)
    #           nm_parts <- character(d)
    #           parts_meta <- vector("list", d)
    #           for (j in seq_len(d)) {
    #             chosen_idx[[j]] <- prim_d[[j]]$idxs[[sel[j]]]
    #             nm_parts[j] <- prim_d[[j]]$names[[sel[j]]]
    #             parts_meta[[j]] <- prim_d[[j]]$meta[[sel[j]]]
    #           }
    #           idx_prod <- interN(chosen_idx)
    #           if (!length(idx_prod)) next
    #           add_cols(
    #             list(idx_prod),
    #             paste0(nm_parts, collapse = ":"),
    #             list(list(type = "interaction", order = d, parts = parts_meta))
    #           )
    #         }
    #       }
    #     }
    #   }
    #
    #   if (chunk_used == 0L) stop("No basis columns created.")
    #
    #   I_chunks <- I_chunks[seq_len(chunk_used)]
    #   J_chunks <- J_chunks[seq_len(chunk_used)]
    #   name_chunks <- name_chunks[seq_len(chunk_used)]
    #   colmeta_chunks <- colmeta_chunks[seq_len(chunk_used)]
    #
    #   i_vec <- unlist(I_chunks, use.names = FALSE)
    #   j_vec <- unlist(J_chunks, use.names = FALSE)
    #   names_all <- unlist(name_chunks, use.names = FALSE)
    #
    #   X <- Matrix::sparseMatrix(
    #     i = if (length(i_vec)) i_vec else integer(),
    #     j = if (length(j_vec)) j_vec else integer(),
    #     x = if (length(i_vec)) 1L else integer(),
    #     dims = c(n, if (length(j_vec)) max(j_vec) else 0L),
    #     dimnames = list(NULL, names_all)
    #   )
    #
    #   list(
    #     X = X,
    #     colnames = names_all,
    #     meta = list(
    #       per_var = meta_per_var,
    #       per_col = colmeta_chunks
    #     )
    #   )
    # },


    update_cross_validation_argument= function(nfold){

      .self$fit_arguments[['nfolds']] <- nfold

    },

    hal_prepare_new = function(DT_new, hal_obj) {
      stopifnot(requireNamespace("data.table", quietly = TRUE))
      stopifnot(requireNamespace("Matrix", quietly = TRUE))

      if (!data.table::is.data.table(DT_new)) DT_new <- data.table::as.data.table(DT_new)

      vars <- names(hal_obj$meta$per_var)
      if (!all(vars %in% names(DT_new))) stop("DT_new is missing required variables")

      DT_local <- DT_new[, ..vars]
      # Coerce char/logical to factor for consistency
      for (v in vars) {
        x <- DT_local[[v]]
        if (is.character(x) || is.logical(x)) DT_local[[v]] <- factor(x)
      }

      n <- nrow(DT_local)
      colnames_tr <- hal_obj$colnames
      colmeta <- hal_obj$meta$per_col
      if (length(colmeta) != length(colnames_tr)) stop("hal_obj meta mismatch")

      # Helper: two-pointer intersection on sorted integer vectors
      inter2 <- function(a, b) {
        i <- 1L; j <- 1L; la <- length(a); lb <- length(b)
        if (!la || !lb) return(integer())
        out <- integer(min(la, lb)); k <- 0L
        while (i <= la && j <= lb) {
          ai <- a[i]; bj <- b[j]
          if (ai < bj) { i <- i + 1L
          } else if (ai > bj) { j <- j + 1L
          } else { k <- k + 1L; out[k] <- ai; i <- i + 1L; j <- j + 1L }
        }
        if (k) out[seq_len(k)] else integer()
      }
      interN <- function(lst) {
        if (length(lst) == 1L) return(lst[[1L]])
        acc <- inter2(lst[[1L]], lst[[2L]])
        if (!length(acc) || length(lst) == 2L) return(acc)
        for (i in 3L:length(lst)) {
          acc <- inter2(acc, lst[[i]])
          if (!length(acc)) break
        }
        acc
      }

      # Precompute, per variable, all primitives needed anywhere in per_col
      need_num <- lapply(vars, function(v) numeric())
      names(need_num) <- vars
      need_fac <- lapply(vars, function(v) character())
      names(need_fac) <- vars

      for (m in colmeta) {
        if (m$type == "main") {
          if (m$kind == "numeric") need_num[[m$var]] <- unique(c(need_num[[m$var]], m$cutpoint))
          else                     need_fac[[m$var]] <- unique(c(need_fac[[m$var]], m$level))
        } else if (m$type == "interaction") {
          for (p in m$parts) {
            if (p$kind == "numeric") need_num[[p$var]] <- unique(c(need_num[[p$var]], p$cutpoint))
            else                     need_fac[[p$var]] <- unique(c(need_fac[[p$var]], p$level))
          }
        } else stop("Unknown column meta type")
      }

      # Build index maps: var -> key -> sorted row indices
      num_map <- vector("list", length(vars)); names(num_map) <- vars
      fac_map <- vector("list", length(vars)); names(fac_map) <- vars

      for (v in vars) {
        x <- DT_local[[v]]
        # numeric
        if (is.numeric(x) && length(need_num[[v]])) {
          cuts <- sort(unique(need_num[[v]]))
          nn <- which(!is.na(x))
          if (length(nn)) {
            ord <- nn[order(x[nn])]
            xs <- x[ord]
            pos <- findInterval(cuts, xs, rightmost.closed = TRUE, all.inside = TRUE)
            idxs <- lapply(pos, function(k) if (k > 0L) ord[seq_len(k)] else integer())
          } else {
            idxs <- lapply(cuts, function(k) integer())
          }
          names(idxs) <- formatC(cuts, digits = 15, format = "fg")  # stable keys
          num_map[[v]] <- idxs
        }
        # factor
        if (is.factor(x) && length(need_fac[[v]])) {
          chr <- as.character(x)
          lv_needed <- need_fac[[v]]
          idxs <- lapply(lv_needed, function(L) which(!is.na(chr) & chr == L))
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
            key <- formatC(m$cutpoint, digits = 15, format = "fg")
            idx <- num_map[[m$var]][[key]]
          } else {
            idx <- fac_map[[m$var]][[m$level]]
          }
        } else { # interaction
          parts_idx <- vector("list", length(m$parts))
          for (k in seq_along(m$parts)) {
            p <- m$parts[[k]]
            if (p$kind == "numeric") {
              key <- formatC(p$cutpoint, digits = 15, format = "fg")
              parts_idx[[k]] <- num_map[[p$var]][[key]]
            } else {
              parts_idx[[k]] <- fac_map[[p$var]][[p$level]]
            }
            if (is.null(parts_idx[[k]]) || !length(parts_idx[[k]])) { parts_idx <- list(integer()); break }
          }
          idx <- if (length(parts_idx) == 1L) parts_idx[[1L]] else interN(parts_idx)
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

      X <- Matrix::sparseMatrix(
        i = if (length(i_vec)) i_vec else integer(),
        j = if (length(j_vec)) j_vec else integer(),
        x = if (length(i_vec)) 1L else integer(),
        dims = c(n, length(colnames_tr)),
        dimnames = list(NULL, colnames_tr)
      )
      X
    },


    fit = function(data, ...) {

      data_copy = data.table::copy(data)

      group_cols <- c(.self$covariates,.self$treatment)[complete.cases(c(.self$covariates,.self$treatment))]

      data_copy <- data_copy[, .(tij = sum(tij), deltaij = sum(deltaij)), by = c(group_cols, "node", "k")]

      data_copy<-data_copy[complete.cases(data_copy),]

      x_pp <- hal_basis(
        vars = c(group_cols,"node"),
        DT = data_copy,
        max_interaction = .self$max_degree,
        knots_per_order = .self$num_knots
      )

      .self$fit_arguments[['y']] <- as.numeric(data_copy[['deltaij']])

      .self$fit_arguments[['offset']] <-  log(data_copy[['tij']])


      if(!.self$penalise_nodes){

        .self$fit_arguments[['penalty.factor']] <- 1- grepl('^I\\(\\s*node\\s*==[^)]*\\)$', x_pp$colnames)

      }

      .self$fit_arguments[['x']] <- x_pp$X

      .self$covariates_attributes_matrix <-x_pp

      if (.self$cross_validation) {

        prefit_args <- .self$fit_arguments

        if(!is.null(.self$maxit_prefit)){prefit_args[['maxit']] <- .self$maxit_prefit}

        cv_fit <- do.call(cv.glmnet, prefit_args)

        if (is.null(cv_fit$glmnet.fit)) return(cv_fit)

        lambda_grid <- cv_fit$lambda

        if (!all(is.na(.self$lambda_grid))) {
          lambda_grid <- c(.self$lambda_grid, lambda_grid)
          lambda_grid <- lambda_grid[complete.cases(lambda_grid)]
          lambda_grid <- sort(unique(lambda_grid), decreasing = TRUE)
        }

        preds <- predict(cv_fit$glmnet.fit,
                         newx = x_pp$X,
                         newoffset = log(data_copy[['tij']]),
                         type = "response",
                         s = lambda_grid)

        if (is.null(dim(preds))) {
          preds <- matrix(preds, ncol = 1L)
        }

        mu <- pmax(preds, 1e-12)
        y_vec <- as.numeric(data_copy[['deltaij']])
        nll <- -colSums(y_vec * log(mu) - mu)
        lambda_opt <- lambda_grid[which.min(nll)]

        glmnet_args <- .self$fit_arguments
        glmnet_args[['lambda']] <- lambda_opt
        glmnet_args[['nfolds']] <- NULL

        fit <- do.call(glmnet, glmnet_args)
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
        if (is.null(fit)) return(TRUE)
        df <- if (inherits(fit, "cv.glmnet")) fit$glmnet.fit$df else fit$df
        if (!is.null(df)) return(all(df == 0))
        beta <- if (inherits(fit, "cv.glmnet")) fit$glmnet.fit$beta else fit$beta
        if (is.null(beta)) return(TRUE)
        all(beta == 0)
      }

      if (is_empty_model(model)) {
        return(rep(NA_real_, nrow(newdata)))
      }

      X_new <- .self$hal_prepare_new(newdata, .self$covariates_attributes_matrix)



      out <- predict(model,
                     ...,
                     newx=X_new,
                     newoffset = log(newdata[['tij']]),
                     type = "response")

      return(out)


    },

    private_predictor = function(model, newdata, ...) {

      is_empty_model <- function(fit) {
        if (is.null(fit)) return(TRUE)
        df <- if (inherits(fit, "cv.glmnet")) fit$glmnet.fit$df else fit$df
        if (!is.null(df)) return(all(df == 0))
        beta <- if (inherits(fit, "cv.glmnet")) fit$glmnet.fit$beta else fit$beta
        if (is.null(beta)) return(TRUE)
        all(beta == 0)
      }

      if (is_empty_model(model)) {
        return(rep(NA_real_, nrow(newdata)))
      }

      X_new <- .self$hal_prepare_new(newdata, .self$covariates_attributes_matrix)


      out <- predict(model,
                     ...,
                     newx=X_new,
                     newoffset = log(1),
                     type = "response")

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
    id="character"
  ),
  methods = list(

    initialize = function(covariates = NULL,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          intercept = TRUE,
                          add_nodes = TRUE,
                          penalise_nodes = FALSE,
                          recycle_information = FALSE,
                          id=NA_character_,
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

    update_cross_validation_argument= function(nfold){



    },

    save_meta_data= function(id, ...){

      .self$id <- id

    },

    datapp = function(train_data = NULL, validation_data = NULL) {
      if (!is.null(train_data)) {
        train.mf <- model.frame(as.formula(.self$formula),
                                rbind(train_data, validation_data),
                                drop.unused.levels = FALSE)
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

    fit = function(data, ...) {

      .extract_symbols <- function(term) {
        # try as a bare expression, then as a RHS of a formula
        out <- tryCatch(
          all.vars(str2lang(term)),
          error = function(e) {
            tryCatch(
              all.vars(stats::terms(stats::as.formula(paste("~", term)))),
              error = function(e2) character(0)
            )
          }
        )
        out
      }


      group_cols <- c(.self$covariates,.self$treatment)[complete.cases(c(.self$covariates,.self$treatment))]


      group_cols <- group_cols[is.character(group_cols) & !is.na(group_cols) & nzchar(group_cols)]
      group_cols <- unlist(lapply(group_cols, .extract_symbols), use.names = FALSE)

      group_cols <- unique(group_cols)

      data <- data[, .(tij = sum(tij), deltaij = sum(deltaij)), by = c(group_cols, "node", "k")]

      data <- data[complete.cases(data), ]
      .self$fit_arguments$formula <- as.formula(.self$formula)
      .self$fit_arguments$data <- data
      .self$fit_arguments$offset <- log(data[['tij']])
      fit <- do.call(.self$learner, .self$fit_arguments)
      return(fit)
    },


    predictor = function(model, newdata, ...) {

      pred <- predict(model,
                      newdata = newdata,
                      type = "response",
                      offset = log(newdata[['tij']]),
                      ...)
      return(pred)
    },

    private_predictor = function(model, newdata, ...) {


      pred <- predict(model,
                      newdata = newdata[node %in% model$xlevels$node,],
                      type = "response",
                      offset = log(1),#log(newdata[['tij']]),
                      ...)


      if(all((levels(newdata$node) %in% model$xlevels$node))){

        return(pred)

      }else{


        #
        newdata[node %in% model$xlevels$node,predictions_model:=pred]

        return(as.array(newdata$predictions_model))




        }
    }
  )
)
