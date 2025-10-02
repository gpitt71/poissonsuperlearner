#' Poisson Super Learner implementation
#'
#' This function implements the Poisson SuperLearner.
#'
#' @param data \code{data.frame}, input data to be pre-processed.
#' @param id \code{character}, identifier column.
#' @param status \code{character}, status column.
#' @param nodes \code{numeric}, time grid to construct the piece-wise constant model.
#'
#' @return superlearner
#'
#' @export
Superlearner <- function(data,
                         id = "id", #
                         stratified_k_fold = FALSE, #
                         start_time = NULL, #
                         end_time = NULL, #
                         status = "status", #
                         event_time = NULL, #
                         learners,
                         number_of_nodes = NULL, #
                         nodes = NULL,
                         min_depth=2,
                         meta_learner_algorithms = c("glm","glmnet"),
                         variable_transformation = NULL,
                         nfold = 3, #
                         ...) {



  # Multiple checks about the input ----

  ############
  check_1 <- is.null(start_time) & !is.null(end_time)
  check_2 <- !is.null(start_time) & is.null(end_time)
  check_3 <- (!is.null(start_time) ||
                !is.null(end_time)) & !is.null(event_time)

  # give names to the learners
  if(is.null(names(learners))){

    names(learners) <- paste0("learner_",1:length(learners))
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

  if (min_depth < 2) {
    warning("The minimum number of learners to build an ensemble is two: min_depth will be set to two.")

  }

  ############


  # Data pre-processing ----

  # save some relevant values
  n <- length(unique(data[[id]]))#nrow(data)




  if (!(0 %in% data[[status]])) {
    warning(
      paste0(
        "There is no value of ",
        status,
        " equal to zero. We will consider the data uncensored."
      )
    )
    n_crisks <- length(unique(data[[status]]))
    uncensored_01 <-TRUE

  } else{
    n_crisks <- length(unique(data[[status]])) - 1
    uncensored_01 <-FALSE
  }

  # Pre-process the data

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
      grid_nodes <- seq(min(data[[event_time]]), max(data[[event_time]]) + 1, length.out = as.integer(number_of_nodes))

    } else{
      if (is.null(nodes)) {
        grid_nodes <- sort(unique(data[[event_time]]))



        # grid_nodes <- grid_nodes[-((length(grid_nodes) - 2):length(grid_nodes))]

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
      uncensored_01=uncensored_01
    )



  }



  lhs_string = NULL


  ## Transform the variables if needed ----
  if (!is.null(variable_transformation)) {


    apply_transformations(dt,variable_transformation)


    # Take the variable that we transform
#
#     lhs_vars <- trimws(unlist(strsplit(
#       strsplit(variable_transformation, "~")[[1]][1], "\\+"
#     )))
#     lhs_string <- paste(lhs_vars, collapse = ", ")
#
#     # Take the transformation
#     rhs_vars <- trimws(unlist(strsplit(
#       strsplit(variable_transformation, "~")[[1]][2], "\\+"
#     )))
#     rhs_string <- paste(rhs_vars, collapse = ", ")
#
#
#     eval(parse(
#       text = paste0("
#                dt[,c('", lhs_string
#
#                     , "'):=list(", rhs_string
#                     , ")]
#                ")
#     ))



  }


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

  ## Exploit Poisson likelihood and try to simplify (where possible) the covariates combinations to make the implementation faster ----

  columns_of_interest <- unlist(lapply(learners, function(x) {
    vars <- unique(c(x$covariates, x$treatment))
    vars <- vars[is.character(vars) & !is.na(vars) & nzchar(vars)]
    unlist(lapply(vars, .extract_symbols), use.names = FALSE)
  }), use.names = FALSE)


  columns_of_interest <- unique(columns_of_interest[(complete.cases(columns_of_interest))])

  dt <- dt[, .(tij = sum(tij), deltaij = sum(deltaij),number_of_observations_tmp=.N,id=max(id)), by = c(unique(c(columns_of_interest, lhs_string)), "node", "k")]

  dt[,number_of_observations_tmp:=NULL]

  #### fino a qui.

  ## Splitting in folds ----

  if (stratified_k_fold) {
    setDT(data)

    dt_id <- eval(parse(text = paste0(
      "data[,last(", status, "),by = ", id, "]"
    )))

    eval(parse(text = paste0(
      "setnames(dt_id,'V1','", status, "')"
    )))

    dt_id <- stratified_sampling(dt_id, id, status, nfold)

    dt_id[order(id)]


  } else{
    id_fold <- sample(1:nfold,
                      n,
                      replace = TRUE,
                      prob = rep(1 / nfold, nfold))

    dt_id <- data.table(folder = id_fold, id = unique(data[[id]]))
  }


  dt <- merge(dt, dt_id, by = "id", all.x = T)


  dt_z <- vector("list", n_crisks)

  z_covariates <- paste0("Z", 1:length(learners))



  # if only one learner is present, we simply perform a CV ----
  if (length(learners) == 1) {
    warning("Only one learner was provided.")

    if (learners[[1]]$cross_validation == TRUE) {

      message("\n Cross-validated for the learner is performed on the data.")

      learners[[1]]$update_cross_validation_argument(nfold)

    } else{
      warning(
        "\n The learner was provided with given hyper-parameters. It will be applied with the given configuration."
      )

    }

    # The learner on the full dataset ----
    training_data <- split(dt, by = "k")
    learner_fit <- mapply(function(x){

      out <- learners[[1]]$fit(x)
      return(out)
    },
    training_data,
    SIMPLIFY = FALSE
      )



    #lapply(learners, function(f) f$fit(dt))

    # The learner on the full dataset ----
    fitted_values <- mapply(function(x,newdata){
      learners[[1]]$predictor(model=x,
                              newdata=newdata)
    },
    learner_fit,
    MoreArgs = list(newdata=dt),
    SIMPLIFY = FALSE
    )



    one_learner_out <- list()

    for(causes in 1:n_crisks){


      one_learner_out[[causes]] <- list(
        model = NULL,
        learners_fit=learner_fit[[causes]],
        meta_learner_fit = NULL,
        fitted_values = fitted_values[[causes]]
      )

    }

    #lapply(learners, function(f) f$predictor(model= model, newdata = newdata))


    out <- list(
      learners = learners,
      metalearner = NULL, # it does not exict in this scenario
      superlearner = one_learner_out,
      data_info = list(
        id = id,
        status = status,
        event_time = event_time,
        start_time = start_time,
        end_time = end_time,
        nodes = sort(unique(as.numeric(levels(dt$node)))),
        nfold = nfold,
        maximum_followup = maximum_followup,
        n_crisks=n_crisks,
        variable_transformation = variable_transformation,
        interval_data_type = interval_data_type
      )
    )

    class(out) <- "poisson_superlearner"

    return(out)

  }



  # Train your models in the training set ----
  for (ix in 1:nfold) {


    # Training data for each competing risk ----
    tmp_train <- dt[folder != ix, ]
    training_data <- split(tmp_train, by = "k")
    # Validation data for each competing risk ----
    tmp_val <- dt[folder == ix, ]
    validation_data <- split(tmp_val, by = "k")

    # dt[, train_01:=folder != ix]
    # fd <- split(dt, by = "k")

    # we find the pseudo observations for each fold ----
    pseudo_observations <- mapply(
      function(training_data,
               validation_data,
               competing_risk,
               learners,
               z_covariates,
               ix)
        create_pseudo_observations(training_data, validation_data, competing_risk, learners, z_covariates, ix),
      # create_pseudo_observations(data, learners, z_covariates, ix),
      # fd,
      training_data = training_data,
      validation_data = validation_data,
      competing_risk=as.list(1:n_crisks),
      MoreArgs = list(
        learners = learners,
        ix = ix,
        z_covariates = z_covariates
      ),
      SIMPLIFY = FALSE
    )



    dt_z <- mapply(function(x, y)
      rbind(x, y), dt_z, pseudo_observations, SIMPLIFY = FALSE)

  }


  data_by_competing_risk <- split(dt, by = "k")

  # Compute oos deviance scores for the learners ----

  dt_learners <- rbindlist(
    lapply(seq_along(dt_z), function(i) {
      DT <- copy(dt_z[[i]])
      # detect all Z* columns (flexible: Z1, Z_2, Zabc3, etc.)
      zcols <- grep("^Z", names(DT), value = TRUE)
      melt(
        DT,
        id.vars = c("id", "folder", "node"),
        measure.vars = zcols,
        variable.name = "learner",
        value.name = "pwch"
      )[, which := paste0("pwch_", i)]
    }),
    use.names = TRUE, fill = TRUE
  )

  # transform back to exponential
  dt_learners[,pwch:=exp(pwch)]

  dt_learners[, learner_idx := fifelse(
    grepl("\\d+", learner),
    as.integer(gsub("\\D+", "", learner)),
    match(learner, unique(learner))  # stable ordering for non-numeric suffixes
  )]
  dt_learners[, learner := paste0("learner_", learner_idx)][, learner_idx := NULL]

  dt_learners <- dcast(
    dt_learners,
    id + folder + node + learner ~ which,
    value.var = "pwch"
    # , fun.aggregate = mean  # uncomment if duplicates exist within a cell
  )

  dt_learners[,node:=factor(node,levels = sort(as.numeric(levels(node))),ordered=TRUE)]

  setorder(dt_learners, id, node, learner)


  # Actual deviance computation

  pwch_cols <- paste0("pwch_",1:n_crisks)

  # save sum of pwch

  sum_of_hazards <- paste(pwch_cols, collapse = " + ")

  pwch_dot_string <- paste0("dt_learners[, pwch_dot :=",sum_of_hazards,"]")

  eval(parse(text = pwch_dot_string))

  dt_learners<- merge(dt_learners, dt[k==1,.(id,node,tij)], by = c("id", "node"))

  # compute cumulative hazard

  mapply(function(pwch, name) {
    dt_learners[, (paste0("cumulative_hazard_", name)) := cumsum(get(pwch) * tij), by = .(id,learner)]
  }, pwch_cols, gsub("pwch_", "", pwch_cols))


  # compute survival function

  hazard_terms <- paste0("cumulative_hazard_", 1:n_crisks)
  sum_expr <- paste(hazard_terms, collapse = " + ")
  survival_function_string <- paste0("dt_learners[, survival_function := exp(-(", sum_expr, "))]")

  eval(parse(text = survival_function_string))

  mapply(function(pwch, name) {
    dt_learners[, (paste0("hazard_times_time_", name)) := (get(pwch) * tij)]
  }, pwch_cols, gsub("pwch_", "", pwch_cols))


  # get the deltas
  delta_list =mapply(function(risk,data){data[k==risk,][,c("id","node",paste0("delta_",risk)):=list(id,node,deltaij)][,.SD,.SDcols=c("id","node",paste0("delta_",risk))]},as.list(1:n_crisks),MoreArgs=list(data=dt),SIMPLIFY = F)

  delta_list <- merge_deltas(delta_list)

  delta_list[,node:=factor(node,levels = sort(as.numeric(levels(node))),ordered=TRUE)]

  dt_learners<- merge(dt_learners,
                      delta_list,
                      by=c("id","node"),
                      all.x = T)


  #first term

  mapply(function(name) {
    dt_learners[,paste0("cr_contribution_",name):=get(paste0("delta_", name))*log(get(paste0("delta_", name))/get(paste0("hazard_times_time_", name)))]
    dt_learners[,paste0("cr_contribution_",name):=fifelse(is.nan(get(paste0("cr_contribution_",name))),0,get(paste0("cr_contribution_",name)))]
    dt_learners[,paste0("cr_contribution_",name):=get(paste0("cr_contribution_",name))-(get(paste0("delta_", name))-get(paste0("hazard_times_time_", name)))]

  }, as.list(1:n_crisks))


  # second term

  crc_terms <- paste0("cr_contribution_", 1:n_crisks)
  sum_expr <- paste(crc_terms, collapse = " + ")
  crc_string <- paste0("dt_learners[, cr_contribution_tot := ",sum_expr,"]")

  eval(parse(text = crc_string))


  dt_learners<-dt_learners[,.(deviance_i=sum(cr_contribution_tot)),by=.(id,learner)]

  dt_learners <- merge(dt_learners,dt_id,by="id")

  setkey(dt_learners,NULL)

  dt_learners<-dt_learners[, .(deviance_v = 2*sum(deviance_i, na.rm = TRUE)), by = .(learner,folder)]

  dt_learners<-dt_learners[,.(deviance=mean(deviance_v)),by =learner]


  z_covariates_list <- select_covariate_path(dt_learners, z_covariates, min_depth = min_depth)




  # We do another round of glmnet (or glm) for combining the predictors ----
  ## In the future we can add options for using any algorithm.
  # if (meta_learner_algorithms == "glmnet") {
  #
    # meta_learner <- Learner_glmnet(
    #   covariates = meta_learner_covariates,
    #   cross_validation = nested_cross_validation_meta_learner,
    #   intercept = add_intercept_metalearner,
    #   add_nodes = add_nodes_metalearner,
    #   penalise_nodes = penalise_nodes_metalearner,
    #   ...
    # )
  # } else{
  #   # meta_learner <- Learner_glm(covariates = z_covariates,
  #   #                             add_nodes = add_nodes_metalearner,
  #   #                             intercept = add_intercept_metalearner)
  #   meta_learner <- Metalearner_glm(
  #     covariates = "node:I(Z1-Z2)", #c(z_covariates,paste0(z_covariates,":node")),
  #     intercept = add_intercept_metalearner,
  #     add_nodes = add_nodes_metalearner
  #
  #   )
  #
  # }






#### THIS if-else cycle will be removed for production version of the package: we will only keep the else part.
  if(length(meta_learner_algorithms)==1&!("glmnet"%in%meta_learner_algorithms)&!("glm"%in%meta_learner_algorithms)){

    warning('Only one meta_learner was supplied. Cross-validation on the pseudo-observations will not be performed.')

    # meta_learner <- meta_learner_algorithms

    meta_learner <- meta_learners_candidates(meta_learner_algorithms,
                                              z_covariates)

    meta_learner<-meta_learner[[1]]
    dt_cv_out <- NULL

  }else{


    meta_learners <- lapply(z_covariates_list, function(x) meta_learners_candidates(meta_learner_algorithms,x))

    meta_learners <- unlist(meta_learners, recursive=FALSE)

    # meta_learners <- meta_learners_candidates(meta_learner_algorithms,z_covariates)

    # A second round of cross-validation

    # id_fold <- sample(1:nfold,
    #                   n,
    #                   replace = TRUE,
    #                   prob = rep(1 / nfold, nfold))
    #
    # dt_id <- data.table(folder = id_fold, id = unique(data[[id]]))
    #
    # for (k in seq_along(dt_z)) {
    #   dt_z[[k]][dt_id, on = "id", folder := i.folder]
    #   data_by_competing_risk[[k]][dt_id, on = "id", folder := i.folder]
    # }


    # Fast mapply used for the cr index. The cycle is on the meta learners and the folders (cannot avoid these two).
    ## output data.table

    dt_cv_out <- NULL


    for(meta_l_ix in seq_along(meta_learners)){


      ## Make predictions in each fold

      tmp_cv <- mapply(
        function(dt,
                 dt_z,
                 cr_ix,
                 nfold,
                 meta_learner
                 )
          meta_learner_cross_validation(dt, dt_z,cr_ix,nfold, meta_learner),
        data_by_competing_risk,
        dt_z,
        as.list(1:n_crisks),
        MoreArgs = list(
          nfold = nfold,
          meta_learner = meta_learners[[meta_l_ix]]

        ),
        SIMPLIFY = FALSE
      )


      # Compute the log-likelihood ----
      ## !!!!!!!!!! This needs to be automated!! Now only for two CR ----
      tmp_cv[[1]][, rn := seq_len(.N), by=.(id, folder,node)]
      tmp_cv[[2]][, rn := seq_len(.N), by=.(id, folder,node)]

      tmp_cv <- merge(tmp_cv[[1]], tmp_cv[[2]], by = c("id", "folder", "node","rn"))[
        , rn := NULL][]
      ## !!!!!!!!!!

      tmp_cv[,node:=factor(node,levels = sort(as.numeric(levels(node))),ordered=TRUE)]

      tmp_cv<- tmp_cv[order(id,node),]

      pwch_cols <- paste0("pwch_",1:n_crisks)

      # save sum of pwch

      sum_of_hazards <- paste(pwch_cols, collapse = " + ")

      pwch_dot_string <- paste0("tmp_cv[, pwch_dot :=",sum_of_hazards,"]")

      eval(parse(text = pwch_dot_string))

      tmp_cv<- merge(tmp_cv, dt[k==1,.(id,node,tij)], by = c("id", "node"))


      # compute cumulative hazard

      mapply(function(pwch, name) {
        tmp_cv[, (paste0("cumulative_hazard_", name)) := cumsum(get(pwch) * tij), by = id]
      }, pwch_cols, gsub("pwch_", "", pwch_cols))

      # compute survival function

      hazard_terms <- paste0("cumulative_hazard_", 1:n_crisks)
      sum_expr <- paste(hazard_terms, collapse = " + ")
      survival_function_string <- paste0("tmp_cv[, survival_function := exp(-(", sum_expr, "))]")

      eval(parse(text = survival_function_string))

      mapply(function(pwch, name) {
        tmp_cv[, (paste0("hazard_times_time_", name)) := (get(pwch) * tij)]
      }, pwch_cols, gsub("pwch_", "", pwch_cols))



      #first term

      mapply(function(name) {
        tmp_cv[,paste0("cr_contribution_",name):=get(paste0("delta_", name))*log(get(paste0("delta_", name))/get(paste0("hazard_times_time_", name)))]
        tmp_cv[,paste0("cr_contribution_",name):=fifelse(is.nan(get(paste0("cr_contribution_",name))),0,get(paste0("cr_contribution_",name)))]
        tmp_cv[,paste0("cr_contribution_",name):=get(paste0("cr_contribution_",name))-(get(paste0("delta_", name))-get(paste0("hazard_times_time_", name)))]

      }, as.list(1:n_crisks))

      # second term

      crc_terms <- paste0("cr_contribution_", 1:n_crisks)
      sum_expr <- paste(crc_terms, collapse = " + ")
      crc_string <- paste0("tmp_cv[, cr_contribution_tot := ",sum_expr,"]")

      eval(parse(text = crc_string))


      tmp_cv<-tmp_cv[,.(deviance_i=sum(cr_contribution_tot)),by=id]


      tmp_cv <- merge(tmp_cv,dt_id,by="id")

      setkey(tmp_cv,NULL)


      tmp_cv<-tmp_cv[, .(deviance_v = 2*sum(deviance_i, na.rm = TRUE)), by = folder][,.(deviance=mean(deviance_v))]

      tmp_cv[['meta_learner']] <- names(meta_learners)[meta_l_ix]

      dt_cv_out <- rbind(dt_cv_out,
                         tmp_cv)


    }


    meta_learner <- dt_cv_out[which.min(deviance)][['meta_learner']]


    meta_learner <- meta_learners[[meta_learner]]


    setnames(dt_cv_out,"meta_learner","model")

    }

  meta_learner_fits <- mapply(
    function(dt,
             dt_z,
             meta_learner,
             learners,
             z_covariates)
      fit_meta_learner(dt, dt_z, meta_learner, learners, z_covariates),
    data_by_competing_risk,
    dt_z,
    MoreArgs = list(
      meta_learner = meta_learner,
      learners = learners,
      z_covariates = z_covariates
    ),
    SIMPLIFY = FALSE


  )



  setnames(dt_learners,"learner","model")

  dt_learners[,covariate:=NULL]

  out <- list(
    learners = learners,
    metalearner = meta_learner,
    superlearner = meta_learner_fits,
    meta_learner_cross_validation=rbind(dt_cv_out,
                                        dt_learners),
    data_info = list(
      id = id,
      status = status,
      event_time = event_time,
      start_time = start_time,
      end_time = end_time,
      nodes = sort(unique(as.numeric(levels(dt$node)))),
      nfold = nfold,
      maximum_followup = maximum_followup,
      n_crisks=n_crisks,
      variable_transformation = variable_transformation,
      interval_data_type = interval_data_type
    )
  )

  class(out) <- "poisson_superlearner"

  return(out)



}
