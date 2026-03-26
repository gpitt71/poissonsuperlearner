#' Fit a single base learner
#'
#' Pre-processes subject-level time-to-event data into a long Poisson format on a
#' piecewise-constant time grid, then fits **one** initialized learner object.
#' For competing risks, a separate model is fit for each event type (cause)
#' using the standard cause-specific Poisson likelihood on the long data.
#'
#' @param data `data.frame`. Subject-level input data (one row per subject).
#' @param learner Reference-class learner object (e.g. from [Learner_glmnet()],
#'   [Learner_hal()] or [Learner_gam()]). Must implement a `$private_fit(dt_long)`
#'   method that fits the learner on long Poisson data for one cause.
#' @param id `character(1)`. Name of the subject identifier column. If not found
#'   in `data`, an `id` column is created automatically.
#' @param stratified_k_fold `logical(1)`. Reserved argument for future fold strategy.
#'   Currently ignored.
#' @param status `character(1)`. Name of the event-status column.
#'   Must be coded with `0` = censoring and `1,2,...,K` for event types (causes).
#'   If there is no `0` in `status`, the data are treated as uncensored.
#' @param event_time `character(1)`. Name of the event/censoring time column.
#'   Must be present in `data`.
#' @param number_of_nodes `numeric(1)` or `NULL`. If not `NULL`, constructs a
#'   quantile-based node grid with `number_of_nodes + 1` cut points (including
#'   endpoints), then adds `0` if missing.
#' @param nodes `numeric` or `NULL`. Explicit time-node grid (cut points). If
#'   supplied, `number_of_nodes` is ignored. `0` is added if missing. Nodes
#'   beyond `max(event_time)` are dropped.
#' @param variable_transformation `list`/`character`/`formula` or `NULL`.
#'   Optional transformations applied to the internally created long Poisson data
#'   before fitting (via `apply_transformations()`).
#' @param ... Additional arguments currently ignored.
#'
#' @return An object of class `base_learner`, i.e. a named `list` with:
#'
#' \describe{
#' \item{model}{The **learner object** that was fit (the input `learner`), stored
#'   for later prediction. This contains the learner specification (e.g.,
#'   covariates, tuning parameters).}
#'
#' \item{learner_fit}{A `list` of fitted model objects, **one per cause**.
#'   Its length equals `data_info$n_crisks`. The list is created by splitting the
#'   internally pre-processed long data by cause indicator `k` and calling
#'   `model$private_fit()` on each split.
#'
#'   \itemize{
#'     \item Names typically correspond to the cause labels `"1"`, `"2"`, ..., `"K"`.
#'     \item Each element is **learner-dependent**: e.g. for `Learner_glmnet` it
#'       may be a `"glmnet"` (often wrapped, e.g. `"fishnet"`) fit; for other
#'       learners it will be whatever `$private_fit()` returns.
#'     \item Each fitted object is trained on long Poisson data representing the
#'       piecewise-constant hazard for that cause across the node intervals.
#'   }}
#'
#' \item{data_info}{A `list` of bookkeeping information needed for prediction and
#'   interpretation:
#'   \describe{
#'     \item{id}{Identifier column name used.}
#'     \item{status}{Status column name used.}
#'     \item{event_time}{Event/censoring time column name used.}
#'     \item{nodes}{Numeric vector of node cut points used for the piecewise grid
#'       (includes `0` and is sorted). These are the interval boundaries used in
#'       the long Poisson representation.}
#'     \item{maximum_followup}{`max(data[[event_time]])`.}
#'     \item{n_crisks}{Number of event types (causes) detected.
#'       If censoring is present (`0` in `status`), then `n_crisks = #unique(status) - 1`;
#'       otherwise `n_crisks = #unique(status)`.}
#'     \item{variable_transformation}{The transformation specification passed in
#'       `variable_transformation` (or `NULL`).}
#'   }}
#' }
#'
#' @examples
#' d <- simulateStenoT1(50, competing_risks = TRUE)
#' lrn <- Learner_glmnet(covariates = c("age", "value_LDL"),
#'                       lambda = 0,
#'                       cross_validation = FALSE)
#' bl <- fit_learner(d,
#'                   learner = lrn,
#'                   id = "id",
#'                   status = "status_cvd",
#'                   event_time = "time_cvd",
#'                   number_of_nodes = 4)
#'
#' @export
fit_learner <- function(data,
                    learner,
                    id = "id",
                    stratified_k_fold = FALSE,
                    status = "status",
                    event_time = NULL,
                    number_of_nodes = NULL,
                    nodes = NULL,
                    variable_transformation = NULL,
                    ...) {
  # Multiple checks about the input ----


  if (!(id %in% names(data))) {
    data[["id"]] <- 1:NROW(data)
    id <- "id"
  }


  maximum_followup = max(data[[event_time]])


  n <- length(unique(data[[id]]))

  if (!(0 %in% data[[status]])) {
    warning(
      paste0(
        "There is no value of ",
        status,
        " equal to zero. We will consider the data uncensored."
      )
    )
    n_crisks <- length(unique(data[[status]]))
    uncensored_01 <- TRUE

  } else{
    n_crisks <- length(unique(data[[status]])) - 1
    uncensored_01 <- FALSE
  }


    #  Handle nodes
    ##Either the nodes are given or we take all of the realised times
    if (!is.null(number_of_nodes)) {

      grid_nodes = quantile(
        data[[event_time]],
        probs = seq(0, 1, length.out = as.integer(number_of_nodes) + 1),
        type = 1,
        names = FALSE
      )



    } else{
      if (is.null(nodes)) {
        grid_nodes <- sort(unique(data[[event_time]]))

      } else{
        grid_nodes <- nodes

      }
    }

    # Add zero if missing
    if (!(0 %in% grid_nodes)) {
      grid_nodes <- c(0, grid_nodes)

    }



    grid_nodes <- grid_nodes[grid_nodes <= max(data[[event_time]])]

    # Actual data pp
    dt <- data_pre_processing(
      data = data,
      id = id,
      status = status,
      nodes = grid_nodes,
      event_time = event_time,
      uncensored_01 = uncensored_01
    )






  lhs_string = NULL

  if (!is.null(variable_transformation)) {
    apply_transformations(dt, variable_transformation)
  }

  training_data <- split(dt, by = "k")

  model_fit <- mapply(function(x) {
    out <- learner$private_fit(x)
    return(out)
  }, training_data, SIMPLIFY = FALSE)


  out <- list(
    model=learner,
    learner_fit=model_fit,
    data_info = list(
      id = id,
      status = status,
      event_time = event_time,
      nodes = sort(unique(as.numeric(levels(dt$node)))),
      maximum_followup = maximum_followup,
      n_crisks=n_crisks,
      variable_transformation = variable_transformation    )
  )

  class(out) <- "base_learner"

  return(out)
}
