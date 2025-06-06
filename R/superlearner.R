#' Poisson Super Learner implementation
#'
#' This function implements the Poisson SuperLearner.
#'
#' @param data \code{data.frame}, input data to be pre-processed.
#' @param id \code{character}, identifier column.
#' @param status \code{character}, status column.
#' @param nodes \code{numeric}, time grid to construct the piece-wise constant model.
#'
#' @return superlearner
#'
#' @export
Superlearner <- function(data,
                         id = "id",
                         stratified_k_fold = FALSE,
                         start_time = NULL,
                         end_time = NULL,
                         status = "status",
                         event_time = NULL,
                         learners,
                         number_of_nodes = NULL,
                         nodes = NULL,
                         meta_learner_algorithm = "glmnet",
                         add_nodes_metalearner = TRUE,
                         add_intercept_metalearner = TRUE,
                         matrix_transformation = FALSE,
                         penalise_nodes_metalearner = TRUE,
                         variable_transformation = NULL,
                         nfold = 3) {
  # Multiple checks about interval data
  # browser()
  check_1 <- is.null(start_time) & !is.null(end_time)
  check_2 <- !is.null(start_time) & is.null(end_time)
  check_3 <- (!is.null(start_time) ||
                !is.null(end_time)) & !is.null(event_time)


  if (check_1 || check_2) {
    stop("For interval data, both start_time and end_time are required.")

  }

  if (check_3) {
    stop("Either provide interval data or censored data")

  }

  if (!is.null(start_time) & !is.null(end_time)) {
    interval_data_type = TRUE

    maximum_followup = max(data[[end_time]])

  } else{
    interval_data_type = FALSE

    maximum_followup = max(data[[event_time]])

  }




  # save some relevant values
  n <- length(unique(data[[id]]))#nrow(data)
  n_crisks <- length(unique(data[[status]])) - 1


  # Pre-process the data

  if (interval_data_type) {
    # Compute here the nodes checking
    if (!is.null(number_of_nodes)) {
      observed_times <- c(data[[start_time]], data[[end_time]])
      grid_nodes <- seq(min(observed_times),
                        max(observed_times) + 1,
                        length.out = as.integer(number_of_nodes))

    } else{
      if (is.null(nodes)) {
        observed_times <- c(data[[start_time]], data[[end_time]])
        grid_nodes <- sort(unique(observed_times))
      } else {
        grid_nodes <- sort(nodes)
      }
    }

    if (!(0 %in% grid_nodes)) {
      grid_nodes <- c(0, grid_nodes)
    }

    # Actual data pp
    dt <- data_pre_processing_interval_data(
      data = data,
      id = id,
      status = status,
      start_time = start_time,
      nodes = grid_nodes,
      end_time = end_time
    )




  } else{
    #  Handle nodes
    ##Either the nodes are given or we take all of the realised times
    if (!is.null(number_of_nodes)) {
      grid_nodes <- seq(min(data[[event_time]]), max(data[[event_time]]) + 1, length.out = as.integer(number_of_nodes))

    } else{
      if (is.null(nodes)) {
        grid_nodes <- sort(unique(data[[event_time]]))


        # browser()
        grid_nodes <- grid_nodes[-((length(grid_nodes) - 2):length(grid_nodes))]

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
      event_time = event_time
    )


    # browser()
  }


  if (matrix_transformation) {
    # browser()

    columns_of_interest <- unlist(lapply(learners, function(x) {
      return(unique(c(x$covariates, x$treatment)))
    }))

    columns_of_interest <- unique(columns_of_interest[(complete.cases(columns_of_interest))])

    # Take the variable that we transform

    lhs_vars <- trimws(unlist(strsplit(
      strsplit(variable_transformation, "~")[[1]][1], "\\+"
    )))
    lhs_string <- paste(lhs_vars, collapse = ", ")

    # Take the transformation
    rhs_vars <- trimws(unlist(strsplit(
      strsplit(variable_transformation, "~")[[1]][2], "\\+"
    )))
    rhs_string <- paste(rhs_vars, collapse = ", ")


    eval(parse(
      text = paste0("
               dt[,c('", lhs_string

                    , "'):=list(", rhs_string
                    , ")]
               ")
    ))


    dt <- dt[, .(tij = sum(tij), deltaij = sum(deltaij)), by = c(unique(c(columns_of_interest, lhs_string)), "node", "k")]


    dt[, c("id") := 1:nrow(dt)]



  }


  if (stratified_k_fold) {
    setDT(data)

    dt_id <- eval(parse(text = paste0(
      "data[,last(", status, "),by = ", id, "]"
    )))

    eval(parse(text = paste0(
      "setnames(dt_id,'V1','", status, "')"
    )))

    dt_id <- stratified_sampling(dt_id, id, status, nfold)

    dt_id[order(id)]


  } else{
    id_fold <- sample(1:nfold,
                      n,
                      replace = TRUE,
                      prob = rep(1 / nfold, nfold))

    dt_id <- data.table(folder = id_fold, id = unique(data[[id]]))
  }




  # browser()
  dt <- merge(dt, dt_id, by = "id", all.x = T)



  dt_z <- vector("list", n_crisks)

  z_covariates <- paste0("Z", 1:length(learners))



  # Train your models in the training set ----
  for (ix in 1:nfold) {
    # Training data for each competing risk ----
    tmp_train <- dt[folder != ix, ]
    training_data <- split(tmp_train, by = "k")
    # Validation data for each competing risk ----
    tmp_val <- dt[folder == ix, ]
    validation_data <- split(tmp_val, by = "k")

    # dt[, train_01:=folder != ix]
    fd <- split(dt, by = "k")
    # browser()
    # we find the pseudo observations for each fold ----
    pseudo_observations <- mapply(
      function(training_data,
               validation_data,
               #data,
               learners,
               z_covariates,
               ix)
        create_pseudo_observations(training_data, validation_data, learners, z_covariates, ix),
      # create_pseudo_observations(data, learners, z_covariates, ix),
      # fd,
      training_data = training_data,
      validation_data = validation_data,
      MoreArgs = list(
        learners = learners,
        ix = ix,
        z_covariates = z_covariates
      ),
      SIMPLIFY = FALSE
    )



    dt_z <- mapply(function(x, y)
      rbind(x, y), dt_z, pseudo_observations, SIMPLIFY = FALSE)

  }


  # browser()
  data_by_competing_risk <- split(dt, by = "k")

  # browser()
  # We do another round of glmnet (or glm) for combining the predictors ----
  ## In the future we can add options for using any algorithm.
  if (meta_learner_algorithm == "glmnet") {
    meta_learner <- Learner_glmnet(
      covariates = z_covariates,
      cross_validation = T,
      intercept = add_intercept_metalearner,
      add_nodes = add_nodes_metalearner,
      penalise_nodes = penalise_nodes_metalearner
    )
  } else{
    meta_learner <- Learner_glm(covariates = z_covariates,
                                add_nodes = add_nodes_metalearner,
                                intercept = add_intercept_metalearner)

  }


  meta_learner_fits <- mapply(
    function(dt,
             dt_z,
             meta_learner,
             learners,
             z_covariates)
      fit_meta_learner(dt, dt_z, meta_learner, learners, z_covariates),
    data_by_competing_risk,
    dt_z,
    MoreArgs = list(
      meta_learner = meta_learner,
      learners = learners,
      z_covariates = z_covariates
    ),
    SIMPLIFY = FALSE


  )


  # browser()
  out <- list(
    learners = learners,
    metalearner = meta_learner,
    superlearner = meta_learner_fits,
    data_info = list(
      id = id,
      status = status,
      event_time = event_time,
      start_time = start_time,
      end_time = end_time,
      nodes = grid_nodes,
      nfold = nfold,
      maximum_followup = maximum_followup,
      matrix_transformation = matrix_transformation,
      variable_transformation = variable_transformation,
      interval_data_type = interval_data_type
    )
  )

  class(out) <- "poisson_superlearner"

  return(out)



}
