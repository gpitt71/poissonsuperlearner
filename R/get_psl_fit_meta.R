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
