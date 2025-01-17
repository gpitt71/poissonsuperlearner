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

