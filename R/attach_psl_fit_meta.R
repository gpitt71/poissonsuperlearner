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
