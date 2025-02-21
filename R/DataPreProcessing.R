#' Individual Data Pre-Processing
#'
#' This function pre-processes the data for the application of our models.
#'
#' @param data \code{data.table}, input data to be pre-processed.
#' @param id \code{character}, identifier column.
#' @param status \code{character}, status column.
#' @param nodes \code{numeric}, time grid to construct the piece-wise constant model.
#'
#' @return \code{data.table} containing the pre-processed data.
#'
#'
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

    # ,tij = create_offset_variable(grid_nodes, time_to_event = ",event_time,")[, 2], deltaij = create_response_variable_c_risks(grid_nodes,time_to_event = ",
    # event_time,
    # "delta = ",
    # status,
    # "event_type = k)), by = tmp]")))

  ## Retrieve covariates

  dt_fit <- merge(dt_fit, dt, by = id, all.x = TRUE)

  # out <- list(data = dt_fit,
  #             id,
  #             status,
  #             covariates,
  #             treatment)
  #
  # class(out) <- "PPData"


  return(dt_fit)

}






#
# if(setting == "survival"){
#
#   if(is.null(nodes)){
#
#     grid_nodes <- c(0,sort(unique(data[,event])))
#
#   }else{
#
#     grid_nodes <- nodes
#
#     if(!(0 %in% grid_nodes)){
#
#       grid_nodes <- c(0,grid_nodes)
#
#     }
#
#   }
#
#   dt_fit <- data %>%
#     group_by(id) %>%
#     reframe(
#       node = eval(parse(
#         text = paste0(
#           "create_offset_variable_survival(nodes=grid_nodes, time_to_event=",
#           event,
#           ",delta=",
#           status,
#           ")[,1]"
#         )
#       )),
#       tij = eval(parse(
#         text = paste0(
#           "create_offset_variable_survival(nodes=grid_nodes, time_to_event=",
#           event,
#           ",delta=",
#           status,
#           ")[,2]"
#         )
#       )),
#       deltaij = eval(parse(
#         text = paste0(
#           "create_response_variable_survival(nodes=grid_nodes, time_to_event=",
#           event,
#           ",delta=",
#           status,
#           ")"
#         )
#       ))
#     )
#
#
#   dt_fit <- merge(dt_fit,data,by="id",all.x=T)
# }
