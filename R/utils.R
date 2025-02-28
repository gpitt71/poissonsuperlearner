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
  dt_fit <- do.call(rbind, replicate(n_crisks, dt, simplify = FALSE))
  ## we create an artificial k index. Table specific.
  dt_fit <- dt_fit[, k := rep(1:n_crisks, each = dim(dt)[1])]


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

  dt_fit <- merge(dt_fit, dt, by = id, all.x = TRUE)

  dt_fit[,c("node",
            "k"):=list(as.factor(node),
                       as.factor(k))]

  return(dt_fit)

}


# Learners ----

datapp_glmnet <- function(data, formula) {
  train.mf  <- model.frame(as.formula(formula), data)
  x  <- model.matrix(attr(train.mf, "terms"), data = data)
  y  <- data[['deltaij']]
  offset <- log(data[['tij']])

  out <- list(x = x, y = y, offset = offset)

  return(out)
}


# Other utils ----
create_formula <- function(covariates=NA_character_,
                           treatment=NA_character_,
                           competing_risks=FALSE){



  if (!any(is.na( covariates))) {
  xs <- paste(covariates, collapse = "+")
  }

  if (!is.na( treatment)) {
    xs <- paste(xs, "+", treatment)
  }

  if (competing_risks) {
    xs <- paste(xs, "+ k")
  }

  out <- paste("deltaij ~", xs, "-1+node+offset(log(tij))", sep =
                 "")

  return(out)




}












