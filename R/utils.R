#' Helper functions
#'
#' This script contains the utils functions that are used in the package.
#'
#' @import data.table

# Data preprocessing ----

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

  out <- c(rep(0,l-1),
           p_holder)

  return(out)
}

create_offset_variable <- function(nodes, delta, time_to_event){

  tmp <- c(nodes[nodes < time_to_event],
           first(nodes[nodes >= time_to_event]))

  # if(delta == 1){
  tmp[length(tmp)] <- time_to_event
  # }

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
                                nodes=NULL
){



  setDT(data)

  #  Handle nodes ----
  ##Either the nodes are given or we take all of the realised times
  if (is.null(nodes)) {
    grid_nodes <- sort(unique(data[[event_time]]))

  } else{
    grid_nodes <- nodes

  }

  # Add zero if missing
  if (!(0 %in% grid_nodes)) {
    grid_nodes <- c(0, grid_nodes)

  }


  # Handle competing risks ----
  ## for each of the competing risks (CR) we need to create a table
  n_crisks <- length(unique(data[[status]])) - 1
  ## the CR tables are stuck on top of each other to allow for possible interactions
  dt_fit <- do.call(rbind, replicate(n_crisks, data, simplify = FALSE))
  ## we create an artificial k index. Table specific.
  dt_fit <- dt_fit[, k := rep(1:n_crisks, each = dim(data)[1])]


  # Data Transformation ----
  tmp <- c(id, "k")

  dt_fit <- eval(parse(text = paste("dt_fit[, .(node = create_offset_variable(grid_nodes, time_to_event = ",
                                    event_time,
                                    ")[, 1]",
                                    ", tij = create_offset_variable(grid_nodes, time_to_event = ",
                                    event_time,
                                    ")[,2]",
                                    ", deltaij = create_response_variable_c_risks(grid_nodes,time_to_event = ",
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
            "k"):=list(as.factor(node),
                       as.factor(k))]

  return(dt_fit)

}


data_pre_processing_interval_data <- function(data, id, status, start_time, end_time, nodes = NULL) {
  setDT(data)

  if (is.null(nodes)) {
    observed_times <- c(data[[start_time]], data[[end_time]])
    grid_nodes <- sort(unique(observed_times))
  } else {
    grid_nodes <- sort(nodes)
  }

  if (!(0 %in% grid_nodes)) {
    grid_nodes <- c(0, grid_nodes)
  }

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

    intervals <- create_offset_variable_interval_data(grid_nodes, st, et)
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

# Learners ----

datapp_glmnet <- function(data, formula) {
  train.mf  <- model.frame(as.formula(formula),
                           data,
                           drop.unused.levels = FALSE)

  # browser()
  x  <- model.matrix(attr(train.mf, "terms"), data = data)
  y  <- data[['deltaij']]
  offset <- log(data[['tij']])

  out <- list(x = x, y = y, offset = offset)

  return(out)
}

# Other utils ----
create_formula <- function(covariates=NA_character_,
                           treatment=NA_character_,
                           competing_risks=FALSE,
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

  out <- paste("deltaij ~", xs, "-1+offset(log(tij))", sep =
                 "")

  return(out)




}


create_pseudo_observations <- function(train_data,
                                       validation_data,
                                       learners,
                                       z_covariates,
                                       ix){
  "
  This function creates the pseudo-observations.

  "

  # browser()
  check <- setdiff(unique(validation_data$node),
                   unique(train_data$node))

  if(length(check)>0){

    template_row <- train_data[.N]

    # Step 3: Duplicate and modify the node column
    dummy_rows <- rbindlist(lapply(check, function(lvl) {
      row <- copy(template_row)
      row[, node := factor(lvl, levels = levels(train_data$node))]
      return(row)
    }))

    dummy_rows[,tij:=1e-4][,deltaij:=0]

    train_data <- rbind(train_data,dummy_rows)

  }

  train_list <- lapply(learners, function(f) f$fit(train_data))

  # Predict on the validation set your pseudo-observations ----
  val_list <- mapply(
    function(f, model, newdata)
      f$predictor(model = model, newdata = newdata),
    learners,
    train_list,
    MoreArgs = list(newdata = validation_data)
  )



  val_list<- apply(as.matrix(val_list),
                   MARGIN = 2,
                   log)

  # Name the columns

  colnames(val_list) <- z_covariates

  # dt_z <- rbind(dt_z,
  dt_z <- data.table(val_list)[, c("id", "folder") := list(validation_data$id, ix)]
  #)

  return(dt_z)

}


fit_meta_learner <- function(dt,
                             dt_z,
                             meta_learner,
                             learners,
                             z_covariates){


  #
  setorder(dt_z, id, "folder")
  setorder(dt, id, "folder")

  dt_z <- cbind(dt_z,dt)

  meta_learner_fit <- meta_learner$fit(dt_z)

  # learners on the full dataset
  full_train_list <- lapply(learners, function(f) f$fit(dt))

  step_0_predictions <- mapply(
    function(f, model, newdata)
      f$predictor(model = model, newdata = newdata),
    learners,
    full_train_list,
    MoreArgs = list(newdata = dt)
  )

  step_0_predictions<- apply(as.matrix(step_0_predictions),
                             MARGIN = 2,
                             log)

  # Name the columns
  colnames(step_0_predictions) <- z_covariates



  step_0_predictions <- cbind(as.data.frame.matrix(step_0_predictions), dt[, c("node", "tij","deltaij")])

  setDT(step_0_predictions)

  fitted_values <- meta_learner$predictor(meta_learner_fit,
                                          step_0_predictions)

  out <- list(
    model = meta_learner,
    learners_fit=full_train_list,
    meta_learner_fit = meta_learner_fit,
    fitted_values = fitted_values
  )



  return(out)


}







