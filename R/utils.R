create_offset_variable_survival <- function(nodes, time_to_event){

  tmp <- c(nodes[nodes < time_to_event],
           first(nodes[nodes >= time_to_event]))

  tmp[length(tmp)] <- time_to_event

  tij <- diff(c(0,tmp))

  return(cbind(tmp,tij))
}

create_response_variable_survival <- function(nodes, time_to_event, delta, event_type){


  l <- sum(nodes <= time_to_event)

  out <- c(rep(0,l-1),
           delta)

  return(out)
}
