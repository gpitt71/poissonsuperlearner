#' Individual Data Pre-Processing
#'
#' This function pre-processes the data for the application of our models.
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
                         stratified_k_fold=FALSE,
                         start_time = NULL,
                         end_time = NULL,
                         status = "status",
                         event_time = NULL,
                         learners,
                         nodes = NULL,
                         meta_learner_algorithm = "glmnet",
                         add_nodes_metalearner=TRUE,
                         add_intercept_metalearner=TRUE,
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

  } else{
    interval_data_type = FALSE

  }




  # save some relevant values
  n <- length(unique(data[[id]]))#nrow(data)
  n_crisks <- length(unique(data[[status]])) - 1


  if (stratified_k_fold) {

    setDT(data)

    dt_id <- eval(parse(text = paste0("data[,last(",status,"),by = ",
                                     id,"]")))

    eval(parse(text = paste0("setnames(dt_id,'V1','",status,"')")))

    dt_id<-stratified_sampling(dt_id,id,status,nfold)

    dt_id[order(id)]


  } else{
    id_fold <- sample(1:nfold,
                      n,
                      replace = TRUE,
                      prob = rep(1 / nfold, nfold))

    dt_id <- data.table(folder = id_fold, id = unique(data[[id]]))
  }





  # Pre-process the data

  if (interval_data_type) {
    dt <- data_pre_processing_interval_data(
      data = data,
      id = id,
      status = status,
      start_time = start_time,
      nodes = nodes,
      end_time = end_time
    )




  } else{
    dt <- data_pre_processing(
      data = data,
      id = id,
      status = status,
      nodes = nodes,
      event_time = event_time
    )

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

    # we find the pseudo observations for each fold ----
    pseudo_observations <- mapply(
      function(training_data,
               validation_data,
               learners,
               z_covariates,
               ix)
        create_pseudo_observations(training_data, validation_data, learners, z_covariates, ix),
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

  data_by_competing_risk <- split(dt, by = "k")

  # browser()
  # We do another round of glmnet (or glm) for combining the predictors ----
  ## In the future we can add options for using any algorithm.
  if (meta_learner_algorithm == "glmnet") {
    meta_learner <- Learner_glmnet(
      covariates = z_covariates,
      cross_validation = TRUE,
      intercept=add_intercept_metalearner,
      add_nodes=add_nodes_metalearner
    )
  } else{
    meta_learner <- Learner_glm(covariates = z_covariates,
                                add_nodes=add_nodes_metalearner)

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


  out <- list(
    learners = learners,
    metalearner = meta_learner,
    superlearner = meta_learner_fits,
    data_info = list(
      id = id,
      status = status,
      event_time = event_time,
      start_time=start_time,
      end_time=end_time,
      nodes = nodes,
      nfold = nfold,
      interval_data_type=interval_data_type
    )
  )

  class(out) <- "poisson_superlearner"

  return(out)



}
