#' Helper functions
#'
#' This script contains the utils functions that are used in the package.
#'
#' @import data.table

# Data preprocessing ----

stratified_sampling <- function(dt,id,status,nfold){

  eval(parse(text=paste("setorder(dt,",status,")")))

  proportions_to_tab <- floor(table(dt[[status]])/nfold)

  missing_to_tot <- table(dt[[status]])-proportions_to_tab*nfold

  replicated_list <- replicate(nfold, proportions_to_tab, simplify = FALSE)

  if(sum(missing_to_tot)>0){replicated_list[[length(replicated_list)]] <- replicated_list[[length(replicated_list)]]+missing_to_tot}

  cols <- c(id,status)
  tmp <- dt[,..cols]

  eval(parse(text=paste("setorder(tmp,",status,")")))

  out <- NULL

  for(v in (1:nfold)){

    stratified_v <- sampling::strata(tmp,c(status),size=replicated_list[[v]],method="srswor")

    stratified_v <-getdata(tmp,stratified_v)

    cond <- tmp[[id]]%in%setdiff(tmp[[id]],
            stratified_v[[id]])

    tmp <- tmp[cond,]




    setDT(stratified_v)

    stratified_v[,folder:=v]

    out <- rbind(out,stratified_v)



  }

  setnames(out,id,"id")
  subset_cols<-c("id","folder")
  out<-out[,..subset_cols]
  return(out)


}


create_offset_variable_survival <- function(nodes, delta, time_to_event){

  tmp <- c(nodes[nodes < time_to_event],
           first(nodes[nodes >= time_to_event]))

  if(delta == 1){tmp[length(tmp)] <- time_to_event}

  tij <- diff(c(tmp))

  grid_nodes <- c(nodes[nodes < time_to_event])

  return(cbind(grid_nodes,tij))
}

create_response_variable_survival <- function(nodes, time_to_event, delta, event_type){


  l <- sum(nodes < time_to_event)

  out <- c(rep(0,l-1),
           delta)

  return(out)
}

preprocess_minmax <- function(varData) {
  X <- as.numeric(varData)
  2 * (X - min(X)) / (max(X) - min(X)) - 1
}


preprocess_catdummy <- function(varData, prefix) {

  X <- as.integer(varData)
  n0 <- length(unique(X))
  n1 <- 2:n0
  addCols <- purrr::map(n1, function(x, y) {as.integer(y == x)}, y = X) %>%
    rlang::set_names(paste0(prefix, n1))
  cbind(data, addCols)
}

poisson_nll <- function(y_true, y_pred, ...) {

  # K        <- backend()
  # browser()
  # y_true   <- K$eval(y_true)
  # y_pred   <- K$eval(y_pred)

  out <- -sum(y_true * log(y_pred[2] * y_pred[1]) - y_pred[2] * y_pred[1])

  return(out)

}


create_response_variable_c_risks <- function(nodes, time_to_event, delta, event_type){
  # browser()
  p_holder <- ifelse(delta == event_type, 1, 0)

  l <- sum(nodes < time_to_event)

  out <- c(rep(0, max(0, l - 1)),
           p_holder)

  return(out)
}

create_offset_variable <- function(nodes, delta, time_to_event){

  # browser()
  tmp <- c(nodes[nodes < time_to_event],
           first(nodes[nodes >= time_to_event]))



  if (all(nodes < time_to_event)) {
    tmp <- c(tmp, time_to_event)
  } else{
    tmp[length(tmp)] <- time_to_event
  }


  tij <- diff(c(tmp))

  grid_nodes <- c(nodes[nodes < time_to_event])

  return(cbind(grid_nodes,tij))
}

create_offset_variable_interval_data <- function(nodes, start_time, end_time) {
  nodes <- sort(unique(nodes))
  nodes_in_range <- nodes[nodes >= start_time & nodes < end_time]

  # Add boundaries if needed
  if (!start_time %in% nodes_in_range) nodes_in_range <- c(start_time, nodes_in_range)
  nodes_after <- nodes[nodes >= end_time]
  first_after <- if (length(nodes_after) > 0) min(nodes_after) else end_time
  tmp <- sort(unique(c(nodes_in_range, end_time, first_after)))

  tmp[tmp > end_time] <- end_time
  tmp <- unique(tmp)

  tij <- diff(tmp)
  grid_nodes <- tmp[-length(tmp)]
  return(cbind(grid_nodes, tij))
}

data_pre_processing <- function(data,
                                id,
                                status,
                                event_time,
                                nodes=NULL,
                                uncensored_01=FALSE
){


  # browser()

  setDT(data)

  # Handle competing risks ----
  ## for each of the competing risks (CR) we need to create a table
  n_crisks <- pmax(length(unique(data[[status]])) - 1+uncensored_01,1)
  ## the CR tables are stuck on top of each other to allow for possible interactions
  dt_fit <- do.call(rbind, replicate(n_crisks, data, simplify = FALSE))
  ## we create an artificial k index. Table specific.
  dt_fit <- dt_fit[, k := rep(1:n_crisks, each = dim(data)[1])]


  # Data Transformation ----
  tmp <- c(id, "k")

  dt_fit <- eval(parse(text = paste("dt_fit[, .(node = create_offset_variable(nodes, time_to_event = ",
                                    event_time,
                                    ")[, 1]",
                                    ", tij = create_offset_variable(nodes, time_to_event = ",
                                    event_time,
                                    ")[,2]",
                                    ", deltaij = create_response_variable_c_risks(nodes,time_to_event = ",
                                    event_time,
                                    ", delta=",
                                    status,
                                    ", event_type = k)",
                                    ")",
                                    ", by = .(",
                                    id,
                                    ", k)",
                                    "]")))

  ## Retrieve covariates

  dt_fit <- merge(dt_fit, data, by = id, all.x = TRUE)

  setnames(dt_fit, c(id),c("id"))

  dt_fit[,c("node",
            "k"):=list(factor(node,levels=as.character(nodes)),
                       as.factor(k))]


  # browser()
  dt_fit[,node:=relevel(node,ref=as.character(last(nodes)))]

  return(dt_fit)

}



data_pre_processing_interval_data <- function(data, id, status, start_time, end_time, nodes = NULL) {
  setDT(data)

  n_crisks <- length(unique(data[[status]])) - 1
  output_list <- list()

  # Identify covariate columns dynamically
  core_cols <- c(id, start_time, end_time, status)
  covariate_names <- setdiff(names(data), core_cols)

  for (i in 1:nrow(data)) {
    row_i <- data[i]
    id_val <- row_i[[id]]
    st <- row_i[[start_time]]
    et <- row_i[[end_time]]
    delta <- row_i[[status]]

    intervals <- create_offset_variable_interval_data(nodes, st, et)
    grid_nodes_i <- intervals[, 1]
    tij_i <- intervals[, 2]

    for (k in 1:n_crisks) {
      deltaij <- as.integer((delta == k) & (grid_nodes_i + tij_i >= et))

      # Construct base interval data
      temp_dt <- data.table(
        id = id_val,
        node = grid_nodes_i,
        tij = tij_i,
        deltaij = deltaij,
        k = k
      )

      # Append covariates dynamically
      for (covar in covariate_names) {
        temp_dt[[covar]] <- row_i[[covar]]
      }

      output_list[[length(output_list) + 1]] <- temp_dt
    }
  }



  dt_fit <- rbindlist(output_list)
  # browser()
  # Encode node as N1, N2, ..., based on ordering
  # unique_nodes <- sort(unique(dt_fit$node))
  # node_labels <- paste0("N", seq_along(unique_nodes))
  # node_map <- setNames(node_labels, unique_nodes)
  # dt_fit[, node := factor(paste0("N", match(node, unique_nodes)))]
  # dt_fit[, k := as.factor(k)]

  dt_fit[, c("node", "k") := list(as.factor(node), as.factor(k))]

  return(dt_fit)
}


## Matrix transformation ----
sl_cut <- function(x, breaks, include.lowest = TRUE, right = TRUE) {

  labels <- paste0(head(breaks, -1), "-", tail(breaks, -1) - if (right) 1 else 0)
  # Return labeled factor
  cut(x, breaks = breaks, labels = labels, include.lowest = include.lowest, right = right)
}


# Learners ----

datapp_glmnet <- function(data, formula) {
  train.mf  <- model.frame(as.formula(formula),
                           data,
                           drop.unused.levels = FALSE)

  x  <- model.matrix(attr(train.mf, "terms"), data = data)
  y  <- data[['deltaij']]
  offset <- log(data[['tij']])

  out <- list(x = x, y = y, offset = offset)

  return(out)
}



create_dicretised_data <- function(data,
                                   id,
                                   start_time = NULL, #
                                   end_time = NULL, #
                                   status, #
                                   event_time = NULL, #
                                   number_of_nodes = NULL, #
                                   nodes = NULL, #
                                   variable_transformation){


  check_1 <- is.null(start_time) & !is.null(end_time)
  check_2 <- !is.null(start_time) & is.null(end_time)
  check_3 <- (!is.null(start_time) ||
                !is.null(end_time)) & !is.null(event_time)


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


  # save some relevant values
  n <- length(unique(data[[id]]))#nrow(data)
  n_crisks <- length(unique(data[[status]])) - 1


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


        # browser()
        grid_nodes <- grid_nodes[-((length(grid_nodes) - 2):length(grid_nodes))]

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
      event_time = event_time
    )


    # browser()
  }


  if (matrix_transformation) {
    # browser()

    columns_of_interest <- unlist(lapply(learners, function(x) {
      return(unique(c(x$covariates, x$treatment)))
    }))

    columns_of_interest <- unique(columns_of_interest[(complete.cases(columns_of_interest))])

    # Take the variable that we transform

    lhs_vars <- trimws(unlist(strsplit(
      strsplit(variable_transformation, "~")[[1]][1], "\\+"
    )))
    lhs_string <- paste(lhs_vars, collapse = ", ")

    # Take the transformation
    rhs_vars <- trimws(unlist(strsplit(
      strsplit(variable_transformation, "~")[[1]][2], "\\+"
    )))
    rhs_string <- paste(rhs_vars, collapse = ", ")


    eval(parse(
      text = paste0("
               dt[,c('", lhs_string

                    , "'):=list(", rhs_string
                    , ")]
               ")
    ))


    dt <- dt[, .(tij = sum(tij), deltaij = sum(deltaij)), by = c(unique(c(columns_of_interest, lhs_string)), "node", "k")]


    dt[, c("id") := 1:nrow(dt)]



  }





  return(dt)



}


# Other utils ----
create_formula <- function(covariates=NA_character_,
                           treatment=NA_character_,
                           competing_risks=FALSE,
                           intercept=FALSE,
                           add_nodes=TRUE){



  if (!any(is.na( covariates))) {
  xs <- paste(covariates, collapse = "+")
  }

  if (!is.na( treatment)) {
    xs <- paste(xs, "+", treatment)
  }

  # if (competing_risks) {
  #   xs <- paste(xs, "+ k")
  # }

  if (add_nodes) {
    xs <- paste(xs, "+ node")
  }

  if(!intercept){
    xs <- paste(xs, "-1")
  }

  out <- paste("deltaij ~", xs, "+offset(log(tij))", sep =
                 "")

  return(out)




}

create_formula_gam <- function(covariates=NA_character_,
                           treatment=NA_character_,
                           competing_risks=FALSE,
                           intercept=FALSE,
                           add_nodes=TRUE){



  if (!any(is.na( covariates))) {
    xs <- paste(covariates, collapse = "+")
  }

  if (!is.na( treatment)) {
    xs <- paste(xs, "+", treatment)
  }

  # if (competing_risks) {
  #   xs <- paste(xs, "+ k")
  # }

  if (add_nodes) {
    xs <- paste(xs, "+ node")
  }

  # if(!intercept){
  #   xs <- paste(xs, "-1")
  # }

  out <- paste("deltaij ~", xs, sep =
                 "")

  return(out)




}

create_formula_hal <- function(covariates=NA_character_,
                           treatment=NA_character_,
                           competing_risks=FALSE,
                           intercept=FALSE,
                           add_nodes=TRUE){



  if (!any(is.na( covariates))) {
    xs <- paste(covariates, collapse = "*")
  }

  if (!is.na( treatment)) {
    xs <- paste(xs, "*", treatment)
  }

  # if (competing_risks) {
  #   xs <- paste(xs, "+ k")
  # }

  if (add_nodes) {
    xs <- paste(c(xs, "node"), collapse = "*")
  }

  if(!intercept){
    xs <- paste(xs, "-1")
  }

  out <- paste("deltaij ~", xs, "+offset(log(tij))", sep =
                 "")

  return(out)




}



create_pseudo_observations <- function(#data,
                                       training_data,
                                       validation_data,
                                       learners,
                                       z_covariates,
                                       ix){
  "
  This function creates the pseudo-observations.

  "


  # browser()
  train_list <- lapply(learners, function(f) f$private_fit(training_data))


  # Predict on the validation set your pseudo-observations ----
  val_list <- mapply(
    function(f, model, newdata)
      f$private_predictor(model = model, newdata = newdata),
    learners,
    train_list,
    MoreArgs = list(newdata = validation_data)
  )

  # browser()

  # val_list<- as.matrix(val_list)


  val_list<- apply(as.matrix(val_list),
                 MARGIN = 2,
                 log)


  # Name the columns

  colnames(val_list) <- z_covariates


  dt_z <-  data.table(val_list)[, c("id", "folder","node") := validation_data[,.(id,folder,node)]] #

  return(dt_z)

}


fit_meta_learner <- function(dt,
                             dt_z,
                             meta_learner,
                             learners,
                             z_covariates){


  # setorder(dt_z, id, "folder")
  # setorder(dt, id, "folder")
  #

  dt_z <- merge(dt_z,dt,by=c("id","folder","node"))

  # browser()

  # dt_z[, virtual_seq := seq_len(.N), by = .(id, folder)]
  # dt[, virtual_seq := seq_len(.N), by = .(id, folder)]
  #
  # dt_z <- merge(dt_z,dt,
  #               by=c("id","folder","virtual_seq"))
  #
  # dt_z[, virtual_seq := NULL]

  meta_learner_fit <- meta_learner$fit(dt_z)

  fitted_values <- meta_learner$predictor(meta_learner_fit,
                                          dt_z)

  # learners on the full dataset ----
  # browser()
  full_train_list <- lapply(learners, function(f) f$fit(dt))

  # step_0_predictions <- mapply(
  #   function(f, model, newdata)
  #     f$predictor(model = model, newdata = newdata),
  #   learners,
  #   full_train_list,
  #   MoreArgs = list(newdata = dt)
  # )
  #
  # step_0_predictions<- apply(as.matrix(step_0_predictions),
  #                            MARGIN = 2,
  #                            log)

  # Name the columns
  # colnames(step_0_predictions) <- z_covariates

  # step_0_predictions <- cbind(as.data.frame.matrix(step_0_predictions), dt[, c("node", "tij","deltaij")])
  #
  # setDT(step_0_predictions)
  #
  # fitted_values <- meta_learner$predictor(meta_learner_fit,
  #                                         step_0_predictions)

  out <- list(
    model = meta_learner,
    learners_fit=full_train_list,
    meta_learner_fit = meta_learner_fit,
    fitted_values = fitted_values
  )



  return(out)


}


# super-learner hyper parameters cross-validation ----


meta_learner_cross_validation <- function(dt,
                             dt_z,
                             cr_ix,
                             nfold,
                             meta_learner){





  # model output ---
  out<-NULL

  dt_z <- merge(dt_z,dt,by=c("id","folder","node"))


  for(v_fold_id in 1:nfold){

  meta_learner_fit <- meta_learner$fit(dt_z[folder!=v_fold_id,])

  oos_data <- copy(dt_z[folder == v_fold_id, ])

  # forced tij to be one
  oos_data[,tij:=1]


  fitted_hazard = data.table(
    oos_data[['id']],
    oos_data[['node']],
    oos_data[['deltaij']],
    v_fold_id,
    as.vector(meta_learner$predictor(meta_learner_fit, oos_data))
  )

  setnames(fitted_hazard,c("id","node",paste0("delta_",cr_ix),"folder",paste0("pwch_", cr_ix)))

  out <- rbind(out,fitted_hazard)

  }

  return(out)


}


evaluate_lkh_cv <- function(data,
                            object){






return(NULL)


}



cv_subject_specific_hazard <- function(cause,
                                       object,
                                       newdata,
                                       ...) {


  setDT(newdata)

  tmp <- copy(newdata)


  cond_zero <- 0 %in% tmp[[object$data_info$event_time]]

  cond_times_larger_than_max <- tmp[[object$data_info$event_time]] > object$data_info$maximum_followup


  # if(all(cond_times_larger_than_max)){
  #
  #   warning(paste0("All the entries in the input times are larger than the maximum follow-up: ",
  #                  as.character(object$data_info$maximum_followup)))
  #   d <- NULL
  #
  # }else
  #
  #   {

  #   eval(parse(
  #     text = paste0(
  #       "
  #   vec_dt <- data.table(
  #
  #   ",
  #       object$data_info$event_time,
  #       " = times[times <= object$data_info$maximum_followup]
  # )
  #   "
  #     )
  #   ))

    # tmp[, dummy := 1]
    # vec_dt[, dummy := 1]

    # Merge on dummy to create Cartesian product
    # data_pp <- merge(tmp, vec_dt, by = "dummy", allow.cartesian = TRUE)[, dummy := NULL]



    # if (is.null(data_pp[[object$data_info$id]])) {
    # data_pp[[object$data_info$id]] <- 1:nrow(data_pp)
    # }

    # no problem writing over id
    # data_pp[[object$data_info$id]] <- 1:nrow(data_pp)


    # if (is.null(data_pp[[object$data_info$status]])) {
    # data_pp[[object$data_info$status]] <- 0
    # }

    data_pp <- data_pre_processing(
      newdata,
      id = object$data_info$id,
      status = object$data_info$status,
      event_time = object$data_info$event_time,
      nodes = object$data_info$nodes
    )




    if(object$data_info$matrix_transformation){


      columns_of_interest <- unlist(lapply(object$learners, function(x){return(unique(c(x$covariates,x$treatment)))}))

      columns_of_interest <- unique(columns_of_interest[(complete.cases(columns_of_interest))])

      # Take the variable that we transform

      lhs_vars <- trimws(unlist(strsplit(strsplit(object$data_info$variable_transformation, "~")[[1]][1], "\\+")))
      lhs_string <- paste(lhs_vars, collapse = ", ")

      # Take the transformation
      rhs_vars <- trimws(unlist(strsplit(strsplit(object$data_info$variable_transformation, "~")[[1]][2], "\\+")))
      rhs_string <- paste(rhs_vars, collapse = ", ")


      eval(parse(text = paste0(
        "
               data_pp[,c('", lhs_string

        , "'):=list(", rhs_string
        , ")]
               "
      )))


      # data_pp <- data_pp[,.(tij = sum(tij),
      #             deltaij=sum(deltaij)), by = c(unique(c(columns_of_interest,lhs_string)),"node","k")]
      #
      #
      # data_pp[,c("id"):=1:nrow(data_pp)]



    }


    # Set covariates for metalearner
    z_covariates <- paste0("Z", 1:length(object$learners))


    # Predict on the validation set your pseudo-observations ----
    # data_pp[,deltatime:=tij][,tij:=1]

    learners_predictions <- mapply(
      function(f, model, newdata)
        f$predictor(model = model, newdata = data_pp),
      object$learners,
      object$superlearner[[cause]]$learners_fit,
      MoreArgs = list(newdata = data_pp)
    )

    # browser()

    pseudo_observations_data <- matrix(apply(as.matrix(learners_predictions, nrow=nrow(newdata), ncol=length(z_covariates)), MARGIN = 2, log),
                                       nrow=nrow(data_pp),
                                       ncol=length(z_covariates))



    # Name the columns

    colnames(pseudo_observations_data) <- z_covariates

    setDT(as.data.frame.matrix(pseudo_observations_data))

    dt_pred <- object$superlearner[[cause]]$model$predictor(object$superlearner[[cause]]$meta_learner_fit,
                                                            newdata =
                                                              cbind(pseudo_observations_data, data_pp))




    data_pp[['pwch']] <- dt_pred


    data_pp <- copy(data_pp)
    # data_pp[,survival_function:=pmin(exp(-cumsum(pwch*deltatime)),1), by=.(id)]


    data_pp <- data_pp[,times:=as.numeric(as.character(node))+tij]


    if(cond_zero){

      data_pp[time==0,
              c('pwch',
                'survival_function'):=list(0,1)]



    }



    columns_ss <- unique(c(colnames(newdata),"times","pwch","deltaij","k"))

    d <- data_pp[,..columns_ss]

  # }


    if(any(cond_times_larger_than_max)){
      # browser()

      d[times > object$data_info$maximum_followup ,c('pwch','survival_function'):=list(NA,NA)]


    }





  return(d)




}

learners_hat <- function(crisk_cause,superlearner,newdata,learners){
  # browser()
  learners_predictions <- mapply(
    function(f, model, newdata)
      f$predictor(model = model, newdata = newdata),
    learners,
    superlearner$learners_fit,
    MoreArgs = list(newdata = newdata)
  )

  return(learners_predictions)



}




