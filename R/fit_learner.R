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
#' Fit a single base learner
#'
#' Pre-processes subject-level time-to-event data into the compressed long
#' Poisson representation used by the current learner interface, then fits one
#' initialized learner object separately for each competing cause.
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

  ## ------------------------------------------------------------
  ## Basic checks and bookkeeping
  ## ------------------------------------------------------------

  data <- data.table::as.data.table(data.table::copy(data))

  if (is.null(event_time) || !(event_time %in% names(data))) {
    stop("'event_time' must name a column in 'data'.", call. = FALSE)
  }

  if (!(status %in% names(data))) {
    stop("'status' must name a column in 'data'.", call. = FALSE)
  }

  if (
    !any(grepl("Learner_", class(learner), fixed = TRUE)) ||
    !is.function(learner$private_fit) ||
    !is.function(learner$private_predictor)
  ) {
    stop(
      "'learner' must be an initialized learner object with private_fit() and private_predictor().",
      call. = FALSE
    )
  }

  if (!(id %in% names(data))) {
    data[, id := seq_len(.N)]
    id <- "id"
  }

  maximum_followup <- max(data[[event_time]])

  if (!(0 %in% data[[status]])) {
    warning(
      paste0(
        "There is no value of ",
        status,
        " equal to zero. We will consider the data uncensored."
      ),
      call. = FALSE
    )

    n_crisks <- length(unique(data[[status]]))
  } else {
    n_crisks <- length(unique(data[[status]])) - 1L
  }


  ## ------------------------------------------------------------
  ## Build node grid
  ## ------------------------------------------------------------

  if (!is.null(number_of_nodes)) {
    grid_nodes <- quantile(
      data[[event_time]],
      probs = seq(0, 1, length.out = as.integer(number_of_nodes) + 1L),
      type = 1,
      names = FALSE
    )
  } else if (is.null(nodes)) {
    grid_nodes <- sort(unique(data[[event_time]]))
  } else {
    grid_nodes <- sort(nodes)
  }

  if (!(0 %in% grid_nodes)) {
    grid_nodes <- c(0, grid_nodes)
  }

  grid_nodes <- sort(unique(grid_nodes))
  grid_nodes <- grid_nodes[grid_nodes <= maximum_followup]


  ## ------------------------------------------------------------
  ## Compressed long-data construction
  ## Same structure expected by Learner_glmnet/Learner_hal/Learner_gam
  ## ------------------------------------------------------------

  grid_dt <- data.table::data.table(node = grid_nodes)
  grid_dt[, prev_node_psl := data.table::shift(node)]
  grid_dt[, (event_time) := node]

  dt <- grid_dt[data, on = event_time, roll = Inf]

  event_time_col <- event_time

  dt[
    ,
    node_start := data.table::fifelse(
      get(event_time_col) == node,
      prev_node_psl,
      node
    )
  ]

  dt <- dt[
    !is.na(node_start),
    c("node", "tij") := list(
      node_start,
      get(event_time_col) - node_start
    )
  ]

  data.table::setnames(dt, old = status, new = "deltaij")

  if (!is.null(variable_transformation)) {
    apply_transformations(dt, variable_transformation)
  }


  ## ------------------------------------------------------------
  ## Fit the same learner once per cause
  ## ------------------------------------------------------------

  causes <- seq_len(n_crisks)

  learner_fit <- learner$private_fit_all_causes(
    data = dt,
    causes = causes,
    grid_nodes = grid_nodes
  )

  ## ------------------------------------------------------------
  ## Output
  ## ------------------------------------------------------------

  out <- list(
    model = learner,
    learner_fit = learner_fit,
    data_info = list(
      id = id,
      status = status,
      event_time = event_time,
      nodes = grid_nodes,
      maximum_followup = maximum_followup,
      n_crisks = n_crisks,
      variable_transformation = variable_transformation
    )
  )

  class(out) <- "base_learner"

  out
}
