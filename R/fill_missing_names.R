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
