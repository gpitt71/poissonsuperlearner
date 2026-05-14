#' Fill missing or empty names in a list-like object
#'
#' Developer note: internal helper used to normalize learner and cause-library
#' names before storing them in fitted objects.
#'
#' Assumes `x` can carry names and that `prefix` is a scalar character prefix.
#' Missing, `NA`, or empty names are replaced by `prefix` plus the element index,
#' and all names are made unique with `make.unique()`.
#'
#' @param x List-like object whose names should be completed.
#' @param prefix Character prefix used for generated names.
#' @returns `x` with completed, unique names.
#' @noRd
fill_missing_names <- function(x, prefix) {
  nm <- names(x)

  if (is.null(nm)) {
    nm <- rep.int("", length(x))
  }

  missing_nm <- is.na(nm) | nm == ""

  nm[missing_nm] <- paste0(prefix, seq_along(x))[missing_nm]

  names(x) <- make.unique(nm, sep = "_")
  x
}

# Data preprocessing ----

#' Create a cause-specific long-format response vector
#'
#' Developer note: internal helper for piecewise-constant hazard preprocessing.
#'
#' Assumes `nodes` are sorted time cut points, `time_to_event` is a single event
#' or censoring time, `delta` is the observed status, and `event_type` is the
#' cause being expanded. The returned vector has zero entries before the terminal
#' interval and a terminal indicator for the requested event type.
#'
#' @param nodes Numeric vector of grid nodes.
#' @param time_to_event Numeric scalar event or censoring time.
#' @param delta Observed event status for one subject.
#' @param event_type Cause index currently being encoded.
#' @returns Numeric vector of interval responses for one subject and cause.
#' @noRd
create_response_variable_c_risks <- function(nodes, time_to_event, delta, event_type){


  p_holder <- ifelse(delta == event_type, 1, 0)

  l <- sum(nodes < time_to_event)

  out <- c(rep(0, max(0, l - 1)),
           p_holder)

  return(out)
}

#' Create interval exposures for one subject
#'
#' Developer note: internal helper used while expanding subject-level data to
#' long Poisson data.
#'
#' Assumes `nodes` are sorted cut points on the follow-up scale and
#' `time_to_event` is a single observed time. The `delta` argument is retained
#' for the historical call signature but is not used. The function truncates the
#' final interval at the observed time and returns interval starts with their
#' exposure widths. It can fail if `nodes` do not contain values compatible with
#' the observed time.
#'
#' @param nodes Numeric vector of grid nodes.
#' @param delta Unused event indicator/status argument.
#' @param time_to_event Numeric scalar event or censoring time.
#' @returns Two-column matrix with `grid_nodes` and `tij` interval lengths.
#' @noRd
create_offset_variable <- function(nodes, delta, time_to_event){

  # if(time_to_event==0){
  #
  #   return(cbind(0,0))
  #
  # }else{
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


  # }
}

#' Expand subject-level data into long Poisson data
#'
#' Developer note: internal preprocessing step shared by base learners and the
#' Super Learner. It mutates `data` by reference after coercing it with
#' `data.table::setDT()`.
#'
#' Assumes `data` contains the supplied identifier, status, and event-time
#' columns, and that `nodes` are sorted cut points covering the needed follow-up
#' range. The output contains one row per subject, cause, and at-risk interval,
#' with interval exposures, cause-specific responses, and factor-coded nodes.
#' Failures generally arise from missing columns, incompatible node values, or
#' invalid event/status encodings.
#'
#' @param data Input subject-level `data.frame` or `data.table`.
#' @param id Name of the subject identifier column.
#' @param status Name of the event-status column.
#' @param event_time Name of the event/censoring time column.
#' @param nodes Numeric vector of grid nodes.
#' @param predictions Historical flag, currently unused.
#' @param uncensored_01 Logical/numeric flag adjusting competing-risk counting.
#' @returns A `data.table` in long Poisson format.
#' @noRd
data_pre_processing <- function(data,
                                id,
                                status,
                                event_time,
                                nodes=NULL,
                                predictions=FALSE,
                                uncensored_01=FALSE
){


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
    dt_fit <- dt_fit[ , {
        tte <- .SD[[1L]]   # time-to-event vector
        del <- .SD[[2L]]   # status/delta vector
        off <- create_offset_variable(nodes, time_to_event = tte)
        .(node_start    = off[, 1L],
          node    = off[, 1L],
          tij     = off[, 2L],
          deltaij = create_response_variable_c_risks(nodes,time_to_event = tte,delta = del,event_type = k)
          )
    }, by = c(id, "k"),
    .SDcols = c(event_time, status)]

    ## Retrieve covariates

  dt_fit <- merge(dt_fit, data, by = id, all.x = TRUE)

  setnames(dt_fit, c(id),c("id"))

  # if(predictions){
  #   maxn <-last(nodes)
  #   lvls <- as.character(sort(unique(nodes)))
  #
  #
  # }else{
  #
  #   maxn <- max(dt_fit$node)
  #   lvls <- as.character(sort(unique(dt_fit$node)))
  # }
  #
  #
  #
  # dt_fit[,c("node",
  #           "k"):=list(factor(node, levels=lvls),
  #                      as.factor(k))]
  #
  # dt_fit[,node:=relevel(node,ref=as.character(maxn))]

  node_values <- sort(unique(as.numeric(nodes)))
  node_labels <- paste0("n", seq_along(node_values))
  names(node_labels) <- as.character(node_values)

  # keep numeric node values only for optional debugging / future use
  # dt_fit[, node_value := node]

  dt_fit[, node := factor(
    node_labels[as.character(node)],
    levels = node_labels
  )]

  dt_fit[, k := as.factor(k)]

  # keep the last grid node as the reference level, consistently in fit and predict
  dt_fit[, node := relevel(node, ref = tail(node_labels, 1L))]

  return(dt_fit)

}


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
# Matrix transformation ----

#' Apply derived-column transformations to long data
#'
#' Developer note: internal helper used before fitting or prediction when the
#' caller supplies `variable_transformation`.
#'
#' Accepts a formula, character string, or list of formulas/strings of the form
#' `lhs ~ rhs`. Multiple left- and right-hand side expressions may be separated
#' by `+` and must have the same length. The function mutates `dt` by reference
#' and returns it invisibly. It stops for unsupported transformation types,
#' malformed strings, unequal LHS/RHS lengths, parse errors, or evaluation
#' errors in the data.table environment.
#'
#' @param dt `data.table` to modify by reference.
#' @param variable_transformation Transformation specification or `NULL`.
#' @returns Invisibly returns the modified `dt`.
#' @noRd
apply_transformations <- function(dt, variable_transformation) {
  stopifnot(data.table::is.data.table(dt))

  if (is.null(variable_transformation)) return(invisible(dt))

  # Normalize to a list of items (strings or formulas)
  items <- variable_transformation
  if (inherits(items, "formula")) items <- list(items)
  else if (is.character(items))    items <- as.list(items)
  else if (is.list(items))         items <- items
  else stop("Unsupported 'variable_transformation' type.")

  # Collect all (lhs name, rhs expression) pairs
  lhs_all  <- character()
  rhs_all  <- vector("list", 0L)

  for (it in items) {
    if (inherits(it, "formula")) {
      # formula: lhs ~ rhs
      lhs_raw <- paste(deparse(it[[2L]]), collapse = "")
      rhs_raw <- paste(deparse(it[[3L]]), collapse = "")
    } else if (is.character(it) && length(it) == 1L) {
      parts <- strsplit(it, "~", fixed = TRUE)[[1L]]
      if (length(parts) != 2L) stop("Each transformation must contain a single '~'.")
      lhs_raw <- parts[1L]
      rhs_raw <- parts[2L]
    } else {
      stop("List elements must be formulas or single strings of the form 'lhs ~ rhs'.")
    }

    lhs_vec <- trimws(strsplit(lhs_raw, "\\+")[[1L]])
    rhs_vec <- trimws(strsplit(rhs_raw, "\\+")[[1L]])

    if (length(lhs_vec) != length(rhs_vec)) {
      stop(
        sprintf("LHS and RHS have different lengths in '%s ~ %s' (%d vs %d).",
                lhs_raw, rhs_raw, length(lhs_vec), length(rhs_vec))
      )
    }

    # Parse each RHS into an expression
    rhs_exprs <- lapply(rhs_vec, function(s) parse(text = s)[[1L]])

    lhs_all <- c(lhs_all, lhs_vec)
    rhs_all <- c(rhs_all, rhs_exprs)
  }

  # Evaluate and assign each transformation.
  # Using one-by-one assignment keeps evaluation within data.table's environment
  # so column names are visible and functions come from the calling env.
  for (i in seq_along(lhs_all)) {
    dt[, (lhs_all[i]) := eval(rhs_all[[i]])]
  }

  invisible(dt)
}

# Other utils ----


#' Build a glmnet Poisson formula
#'
#' Developer note: internal formula constructor for the `Learner_glmnet`
#' design matrix path.
#'
#' Assumes covariate names are syntactically valid formula terms. The node term
#' is always added and the exposure offset is included. Invalid names or missing
#' required columns fail downstream during model-matrix construction.
#'
#' @param covariates Character vector of covariate names or `NA_character_`.
#' @returns Character scalar formula for `deltaij`.
#' @noRd
create_formula_glmnet <- function(covariates=NA_character_){


  xs<- NULL

  if (!any(is.na( covariates))) {
    xs <- paste(covariates, collapse = "+")
  }

  xs <- paste(xs, "+ node")


  out <- paste("deltaij ~", xs, "+offset(log(tij))", sep =
                 "")

  return(out)




}

#' Build a GAM formula
#'
#' Developer note: internal formula constructor for `Learner_gam`.
#'
#' Assumes covariate names are valid formula terms. The node term is always
#' added, and no exposure offset is included here because the backend handles it
#' separately. The `competing_risks` argument is retained for compatibility.
#'
#' @param covariates Character vector of covariate names or `NA_character_`.
#' @param competing_risks Historical flag, currently unused.
#' @returns Character scalar formula for `deltaij`.
#' @noRd
create_formula_gam <- function(covariates=NA_character_,
                           competing_risks=FALSE){



  if (!any(is.na( covariates))) {
    xs <- paste(covariates, collapse = "+")
  }

  xs <- paste(xs, "+ node")


  out <- paste("deltaij ~", xs, sep =
                 "")

  return(out)




}

#' Build a HAL interaction formula
#'
#' Developer note: internal formula constructor for HAL-style basis generation.
#'
#' Assumes covariate names are valid formula terms. Covariates and `node` are
#' joined with `*` so interaction terms can be expanded by downstream tooling.
#' The `competing_risks` argument is retained for compatibility. Invalid names
#' fail later during formula parsing or basis construction.
#'
#' @param covariates Character vector of covariate names or `NA_character_`.
#' @param competing_risks Historical flag, currently unused.
#' @param intercept Logical; whether to keep an intercept.
#' @returns Character scalar formula with exposure offset.
#' @noRd
create_formula_hal <- function(covariates=NA_character_,
                           competing_risks=FALSE,
                           intercept=FALSE){



  if (!any(is.na( covariates))) {
    xs <- paste(covariates, collapse = "*")
  }

  # if (!is.na( treatment)) {
  #   xs <- paste(xs, "*", treatment)
  # }

  # if (competing_risks) {
  #   xs <- paste(xs, "+ k")
  # }

  xs <- paste(c(xs, "node"), collapse = "*")


  if(!intercept){
    xs <- paste(xs, "-1")
  }

  out <- paste("deltaij ~", xs, "+offset(log(tij))", sep =
                 "")

  return(out)




}


# Make fit cheaper

#' Create a lightweight failed-fit sentinel
#'
#' Developer note: internal helper used to avoid retaining heavy failed model
#' objects in learner libraries.
#'
#' The returned object records a failure reason and has classes used by
#' `is_failed_fit()`. It has no side effects.
#'
#' @param reason Character description of the failure.
#' @returns A `psl_failed_fit` sentinel object.
#' @noRd
make_failed_fit <- function(reason = NA_character_) {
  structure(
    list(reason = reason),
    class = c("psl_failed_fit", "psl_fit_stub")
  )
}

#' Test whether an object is a failed-fit sentinel
#'
#' Developer note: internal predicate paired with `make_failed_fit()`.
#'
#' @param model Object to test.
#' @returns Logical scalar.
#' @noRd
is_failed_fit <- function(model) {
  inherits(model, "psl_failed_fit")
}

#' Attach package fit metadata to a model object
#'
#' Developer note: internal helper for retaining lightweight bookkeeping on
#' fitted learner objects.
#'
#' Returns `NULL` and failed-fit sentinels unchanged. Otherwise it sets the
#' `psl_meta` attribute on `model`, which may copy the object depending on R's
#' usual object semantics.
#'
#' @param model Fitted model object, `NULL`, or failed-fit sentinel.
#' @param meta Metadata object to store as an attribute.
#' @returns `model` with `psl_meta` attached, or the unchanged input.
#' @noRd
attach_psl_fit_meta <- function(model, meta = NULL) {
  if (is.null(model) || is_failed_fit(model)) {
    return(model)
  }
  attr(model, "psl_meta") <- meta
  model
}

#' Read package fit metadata from a model object
#'
#' Developer note: internal accessor for metadata stored by
#' `attach_psl_fit_meta()`.
#'
#' @param model Fitted model object.
#' @returns The exact `psl_meta` attribute, or `NULL` if absent.
#' @noRd
get_psl_fit_meta <- function(model) {
  attr(model, "psl_meta", exact = TRUE)
}

# Hal helpers ----

#' Build main-effect HAL primitive bases for one variable
#'
#' Developer note: internal HAL helper that delegates index construction to C++
#' routines and records metadata needed for later basis expansion.
#'
#' Numeric variables produce threshold indicators. Factor variables drop the
#' first level as a reference and produce non-reference level indicators.
#' Unsupported variable types stop with an error.
#'
#' @param x Variable vector.
#' @param v Character variable name used in generated basis labels.
#' @param K Integer/numeric limit passed to the numeric C++ basis constructor.
#' @returns List with index sets, names, primitive metadata, and variable metadata.
#' @noRd
mk_main <- function(x, v, K) {
  if (is.numeric(x)) {
    out <- mk_main_numeric_cpp(x, as.integer(K))
    cps <- out$cutpoints
    idxs <- out$idxs

    cps <- as.numeric(cps)
    if (!length(cps)) {
      return(list(idxs = list(), names = character(),
                  prim_meta = list(), var_meta = list(cutpoints = numeric())))
    }

    # nms <- sprintf("I(%s<=%.10g)", v, cps)
    nms <- sprintf("I(%s>=%.10g)", v, cps)

    prim_meta <- lapply(cps, function(cu) list(var = v, kind = "numeric", cutpoint = cu))
    list(
      idxs = idxs,
      names = nms,
      prim_meta = prim_meta,
      var_meta = list(cutpoints = cps)
    )

  } else if (is.factor(x)) {
    lv <- levels(x)
    if (!length(lv)) {
      return(list(idxs = list(), names = character(),
                  prim_meta = list(), var_meta = list(levels = character(), ref = NA_character_)))
    }

    # Drop a reference level (default: first level)
    ref <- lv[1]
    keep <- setdiff(lv, ref)

    # If only one level, nothing to add
    if (!length(keep)) {
      return(list(
        idxs = list(),
        names = character(),
        prim_meta = list(),
        var_meta = list(levels = lv, ref = ref)
      ))
    }

    # Compute idxs for all levels, then drop the reference entry
    out_all <- mk_main_factor_cpp(as.integer(x), length(lv))
    idxs_all <- out_all$idxs

    # Map level -> position (mk_main_factor_cpp is assumed to return in level order)
    pos_keep <- match(keep, lv)

    idxs <- idxs_all[pos_keep]
    nms <- sprintf("I(%s==%s)", v, make.names(keep, unique = TRUE))
    prim_meta <- lapply(keep, function(L) list(var = v, kind = "factor", level = L, ref = ref))

    list(
      idxs = idxs,
      names = nms,
      prim_meta = prim_meta,
      var_meta = list(levels = lv, ref = ref)
    )

  } else {
    stop(sprintf("Unsupported type for %s", v))
  }
}

#' Append sparse HAL basis columns to a construction state
#'
#' Developer note: internal HAL helper for accumulating sparse matrix chunks.
#'
#' Assumes `state` has fields `p`, `chunk_used`, `I_chunks`, `J_chunks`,
#' `name_chunks`, and `colmeta_chunks`. The C++ helper removes empty/duplicate
#' columns and reports which names and metadata to keep. Returns `state`
#' unchanged when no columns are added.
#'
#' @param state Mutable construction-state list.
#' @param idxs_list List of row-index vectors for candidate columns.
#' @param nm_vec Character vector of candidate column names.
#' @param col_meta_list List of column metadata objects.
#' @returns Updated construction-state list.
#' @noRd
add_cols <- function(state, idxs_list, nm_vec, col_meta_list) {
  out <- add_cols_cpp(idxs_list, p_start = state$p)
  ncol_added <- out$ncol
  if (!ncol_added) return(state)

  keep <- out$keep
  nm_vec <- nm_vec[keep]
  col_meta_list <- col_meta_list[keep]

  state$chunk_used <- state$chunk_used + 1L
  state$I_chunks[[state$chunk_used]] <- out$i
  state$J_chunks[[state$chunk_used]] <- out$j
  state$name_chunks[[state$chunk_used]] <- nm_vec
  state$colmeta_chunks[[state$chunk_used]] <- col_meta_list

  state$p <- state$p + ncol_added
  state
}
