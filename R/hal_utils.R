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
