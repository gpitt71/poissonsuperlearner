#' Fit a single base learner
#'
#' Pre-processes time-to-event data into long Poisson format and fits one
#' initialized learner reference-class object.
#'
#' @param data `data.frame`. Input subject-level data.
#' @param learner Reference-class learner object (e.g. from [Learner_glmnet()],
#'   [Learner_hal()] or [Learner_gam()]).
#' @param id `character(1)`. Identifier column name.
#' @param stratified_k_fold `logical(1)`. Reserved argument for fold strategy.
#' @param status `character(1)`. Event-status column name.
#' @param event_time `character(1)`. Event/censoring time column name.
#' @param number_of_nodes `numeric(1)` or `NULL`. Number of quantile-based nodes.
#' @param nodes `numeric` or `NULL`. Explicit time-node grid.
#' @param variable_transformation `list`/`character`/`formula` or `NULL`.
#' @param ... Additional arguments currently ignored.
#'
#' @return Object of class `base_learner` with components `learner`,
#'   `learner_fit`, and `data_info`.
#'
#' @examples
#' d <- simulateStenoT1(150, competing_risks = TRUE)
#' lrn <- Learner_glmnet$new(covariates = c("age", "value_LDL"))
#' bl <- fit_learner(d, learner = lrn, id = "id", status = "status_cvd", event_time = "time_cvd")
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
    id <<- "id"
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
