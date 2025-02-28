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
