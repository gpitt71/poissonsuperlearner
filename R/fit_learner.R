#' Fit Learner to the data
#'
#' This class is used for fitting individual learners to the data.
#'
#' @param data \code{data.frame}, input data to be pre-processed.
#' @param id \code{character}, identifier column.
#' @param stratified_k_fold \code{logical}, if TRUE we stratify the V-folds in order to obtain the same number of competing events in each folder.
#' @param start_time \code{character}, for left-truncated and right-censored data, the starting point of each observation.
#' @param end_time \code{character}, for left-truncated and right-censored data, the end point of each observation.
#' @param status \code{character}, status column.
#' @param event_time \code{character}, time-to-event column in competing risks and survival applications.
#' @param learners \code{list}, list of learners to include in the ensemble. If only one learner is provided, the learner is applied to the full data.
#' @param min_depth \code{numeric}, minimum number of learners included in the ensemble.
#' @param meta_learner_algorithms \code{character}, the ensemble is estimated using one of the given algorithms. If more than one algorithm is provided, the best performing algorithm is selected based on the Poisson Deviance. Currently, options are \code{"glm"} and \code{"glmnet"}.
#' @param variable_transformation \code{list}, variable transformation(s) to apply to the data. Each element of the list is a character \code{"new_variable ~ some_function(variable_in_the_data)"}
#' @param nfold \code{numeric}, number of V-folds to construct the ensemble.
#' @param number_of_nodes \code{numeric}, number of time points sampled from the observed time to construct the nodes. Alternative to \code{nodes}, if both NULL we take all the observed time points as nodes.
#' @param nodes \code{numeric}, time grid to construct the piece-wise constant model. Alternative to \code{number_of_nodes}, if both NULL we take all the observed time points as nodes.
#'
#' @return A \code{poisson_superlearner} object contains the following output
#' \itemize{
#' \item{\code{learners}: \code{list} containing the learners.}
#' \item{\code{metalearner}: \code{numeric} Training negative log likelihood.}
#' \item{\code{superlearner}:  \code{numeric} Validation  negative log likelihood. Not available for COX.}
#' \item{\code{meta_learner_cross_validation}:  \code{data.table} containing the average Poisson Deviance evaluated on the V-Folds for the learners and the meta-learner.}
#' \item{\code{data_info}: \code{list} containing the following meta-data:}
#'    \itemize{
#'    \item{\code{id}: \code{character}, name of the covariate in the data that contains the id variable.}
#'    \item{\code{status}: \code{character}, name of the covariate in the data that contains the status variable.}
#'    \item{\code{event_time}: \code{character}, for competing risks or survival data it contains the name of the covariate containing the event time.}
#'    \item{\code{start_time}: \code{character}, for interval data it contains the name of the covariate containing the starting time.}
#'    \item{\code{end_time}: \code{character}, for interval data it contains the name of the covariate containing the ending time.}
#'    \item{\code{nodes}: \code{numeric}, it contains the time nodes.}
#'    \item{\code{nfold}: \code{integer} denoting the number of V-Folds.}
#'    \item{\code{maximum_followup}: \code{numeric} denoting the maximum follow-up time observed.}
#'    \item{\code{n_crisks}, \code{integer} denoting the number of competing risks.}
#'    \item{\code{variable_transformation}: \code{list} or \code{array}, containing the data variable transformations implemented. }
#'    \item{\code{interval_data_type}: \code{logical}, \code{TRUE} for interval data. }
#'    }
#' \item{\code{hazard_model}: \code{string} chosen hazard model (COX, NN or XGB)}
#' \item{\code{IndividualDataPP}: starting \code{IndividualDataPP} object.}
#' }
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
