#' Token-only HAL Poisson learner (experimental)
#' @export Learner_hal_3
#' @exportClass Learner_hal_3
Learner_hal_3 <- methods::setRefClass(
  "Learner_hal_3", contains = "Learner_hal_2", methods = list(
    initialize = function(covariates = NA_character_, intercept = TRUE,
                          cross_validation = FALSE, max_degree = 2L,
                          num_knots = c(10L,5L), maxit_prefit = NA_real_,
                          standardize = TRUE, ...) {
      if (cross_validation) stop("Learner_hal_3 prototype requires cross_validation=FALSE.")
      if (length(max_degree) != 1L || is.na(max_degree) || max_degree < 1L)
        stop("`max_degree` must be a positive integer.")
      if (length(standardize) != 1L || is.na(standardize))
        stop("`standardize` must be TRUE or FALSE.")
      callSuper(covariates=covariates,intercept=intercept,
                cross_validation=FALSE,max_degree=max_degree,
                num_knots=num_knots,maxit_prefit=maxit_prefit,...)
      # Learner_hal uses glmnet's standardize=TRUE default.  Keep the same
      # default objective, while allowing an explicit FALSE for comparisons.
      .self$fit_arguments$standardize <- isTRUE(standardize)
    },
    private_fit = function(data,cause,grid_nodes,...) {
      .self$private_fit_all_causes(data,causes=cause,grid_nodes=grid_nodes,...)[[1L]]
    },
    private_fit_all_causes = function(data,causes,grid_nodes,...) {
      tryCatch({
        prep <- hal2_prepare_training_data(data,.self$basic_covariates,causes,grid_nodes)
        tokens <- hal3_tokens(prep$data,c(.self$basic_covariates,"node"),
                              .self$max_degree,.self$num_knots)
        pf <- tokens$default_penalty_factor
        args <- .self$fit_arguments
        if (is.null(args$lambda)) stop("Supply an explicit decreasing lambda sequence.")
        supported <- c("lambda","tol","maxit_irls","maxit_cd","active_set",
                       "strong_rules","rescale_penalty_factor","standardize")
        unknown <- setdiff(names(args),c(supported,"intercept"))
        if(length(unknown))stop("Unsupported options: ",paste(unknown,collapse=", "))
        ans <- setNames(vector("list",length(causes)),as.character(causes))
        for(i in seq_along(causes)) {
          ans[[i]] <- do.call(hal3_fit_tokens,c(list(
            tokens=tokens,y=prep$y[,i],offset=log(prep$data$tij),
            penalty_factor=pf,intercept=.self$intercept),
            args[intersect(names(args),supported)]))
        }
        ans
      },error=function(e) setNames(rep(list(make_failed_fit(conditionMessage(e))),
                                      length(causes)),as.character(causes)))
    },
    private_predictor = function(model,newdata,grid_nodes,...) {
      ix <- match(newdata$node,grid_nodes)
      if(anyNA(ix))stop("Unexpected time node")
      if (is_failed_fit(model)) return(rep(NA_real_,sum(ix)))
      dots <- list(...)
      selected_lambda <- dots[["s"]]
      if (is.null(selected_lambda)) selected_lambda <- model$selected_lambda
      if (length(selected_lambda) != 1L || !is.finite(selected_lambda))
        stop("`s` must select one finite lambda value.")
      expanded <- newdata[rep(seq_len(nrow(newdata)),times=ix)]
      expanded$node <- factor(sequence(ix),levels=seq_along(grid_nodes),
                              labels=paste0("n",seq_along(grid_nodes)))
      hal3_predict_tokens(model,expanded,s=selected_lambda,
                          type="response",newoffset=rep(0,nrow(expanded)))
    }
  )
)
