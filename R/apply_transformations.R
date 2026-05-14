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
