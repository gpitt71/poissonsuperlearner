#' Expand validation rows to a prediction skeleton
#'
#' Developer note: internal helper for rebuilding all intervals up to each
#' validation subject's terminal interval.
#'
#' Assumes `valid_data` is a `data.table`, `grid_nodes` are sorted numeric nodes,
#' and `valid_data[[node_col]]` matches values in `grid_nodes` exactly. The
#' function copies the requested columns, adds prediction row identifiers, and
#' creates interval-level `node`, `tij`, and `deltaij` columns. It stops if a
#' terminal node cannot be matched to `grid_nodes`.
#'
#' @param valid_data Validation data in terminal long-data form.
#' @param grid_nodes Sorted numeric grid nodes.
#' @param cause Cause index used to mark terminal events.
#' @param node_col Name of the terminal node column.
#' @param tij_col Name of the terminal exposure column.
#' @param event_col Name of the terminal event/status column.
#' @param keep_cols Additional columns to retain.
#' @returns A copied `data.table` expanded to all intervals up to terminal time.
#' @noRd
make_validation_skeleton <- function(valid_data,
                                     grid_nodes,
                                     cause,
                                     node_col  = "node",
                                     tij_col   = "tij",
                                     event_col = "deltaij",
                                     keep_cols = NULL) {
  stopifnot(data.table::is.data.table(valid_data))
  stopifnot(is.numeric(grid_nodes))
  stopifnot(!is.unsorted(grid_nodes))

  cols <- unique(c(keep_cols, node_col, tij_col, event_col))

  x <- data.table::copy(valid_data[, .SD, .SDcols = cols])

  data.table::setnames(x, node_col,  "terminal_node")
  data.table::setnames(x, tij_col,   "terminal_tij")
  data.table::setnames(x, event_col, "terminal_event")

  x[, pred_id := .I]

  x[, terminal_node_ix := match(terminal_node, grid_nodes)]

  if (anyNA(x$terminal_node_ix)) {
    stop("Some terminal_node values are not exactly in grid_nodes.")
  }

  widths <- c(diff(grid_nodes), NA_real_)

  idx <- rep.int(seq_len(nrow(x)), times = x$terminal_node_ix)

  out <- x[idx]

  out[, node_ix := sequence(x$terminal_node_ix)]
  out[, node := grid_nodes[node_ix]]

  # Full-width exposure before the terminal node.
  out[, tij := widths[node_ix]]

  # Observed terminal exposure at the terminal node.
  out[node_ix == terminal_node_ix, tij := terminal_tij]

  # Events only occur at the terminal node, and only if event cause == cause.
  out[, deltaij := 0.0]
  out[node_ix == terminal_node_ix, deltaij := as.numeric(terminal_event == cause)]

  out[, `:=`(
    terminal_tij = NULL,
    terminal_event = NULL
  )]

  out[]
}
