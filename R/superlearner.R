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
                         id="id",
                         start_time=NULL,
                         end_time=NULL,
                         status="status",
                         event_time="time",
                         learners,
                         nodes=NULL,
                         nfold = 3){

  # Multiple checks about interval data

  check_1 <- is.null(start_time) & !is.null(end_time)
  check_2 <- !is.null(start_time) & is.null(end_time)
  check_3 <- (!is.null(start_time) || !is.null(end_time)) & !is.null(event_time)


  if(check_1 || check_2){

    stop("For interval data, both start_time and end_time are required.")

  }

  if(check_3){

    stop("Either provide interval data or censored data")

  }

  if(!is.null(start_time) & !is.null(end_time)){

    interval_data_type = TRUE

  }



  # save some relevant values
  n <- length(unique(data[[id]]))#nrow(data)
  n_crisks <- length(unique(data[[status]])) - 1

  id_fold <- sample(1:nfold,n,replace=TRUE,prob = rep(1/nfold,nfold))

  dt_id <- data.table(
    folder = id_fold,
    id = unique(data[[id]])
  )



  # Pre-process the data

  if(interval_data_type){

    dt <- data_pre_processing_interval_data(
      data = data,
      id = id,
      status = status,
      start_time = start_time,
      end_time = end_time
    )


  }else{

    dt <- data_pre_processing(
      data = data,
      id = id,
      status = status,
      nodes = nodes,
      event_time = event_time
    )

  }



  dt <- merge(dt,
              dt_id,
              by="id",
              all.x = T)




  # browser()
  dt_z <- vector("list",n_crisks)

  z_covariates <- paste0("Z", 1:length(learners))
  # browser()
  # Train your models in the training set ----
  for(ix in 1:nfold){

    # browser()

    # Training data for each competing risk ----
    tmp_train <- dt[folder != ix,]
    training_data <- split(tmp_train, by="k")
    # Validation data for each competing risk ----
    tmp_val <- dt[folder == ix,]
    validation_data <- split(tmp_val, by="k")

    # we find the pseudo observations for each fold ----
    pseudo_observations <- mapply(
      function(training_data,
               validation_data,
               learners,
               z_covariates,
               ix)
        create_pseudo_observations(training_data,
                                 validation_data,
                                 learners,
                                 z_covariates,
                                 ix),
      training_data=training_data,
      validation_data = validation_data,
      MoreArgs = list(learners = learners,
                      ix=ix,
                      z_covariates=z_covariates),
     SIMPLIFY = FALSE
    )

    dt_z <- mapply(function(x,y) rbind(x,y),
                   dt_z,
      pseudo_observations,
      SIMPLIFY = FALSE)


    # train_list <- lapply(learners, function(f) f$fit(tmp_train))

  # Predict on the validation set your pseudo-observations ----
    # val_list <- mapply(
    #   function(f, model, newdata)
    #     f$predictor(model = model, newdata = newdata),
    #   learners,
    #   train_list,
    #   MoreArgs = list(newdata = tmp_val)
    # )



    # val_list<- apply(as.matrix(val_list),
    # MARGIN = 2,
    # log)

    # Name the columns

    # colnames(val_list) <- z_covariates
    #
    # dt_z <- rbind(dt_z,
    #               data.table(val_list)[,c("id",
    #                                       "folder"):=list(tmp_val$id,
    #                                                       ix)])


  }

  # browser()
  data_by_competing_risk <- split(dt, by="k")

  # We do another round of glmnet for combining the predictors ----
  meta_learner <- Learner_glmnet(covariates = z_covariates,
                          cross_validation=TRUE,
                          add_nodes=FALSE)


  meta_learner_fits <- mapply(
    function(dt,
             dt_z,
             meta_learner,
             learners,
             z_covariates)
      fit_meta_learner(dt,
                       dt_z,
                       meta_learner,
                       learners,
                       z_covariates),
    data_by_competing_risk,
    dt_z,
    MoreArgs = list(meta_learner=meta_learner,
                    learners=learners,
                    z_covariates=z_covariates),
    SIMPLIFY = FALSE


  )


  # browser()

  # setorder(dt_z, id, "folder")
  # setorder(dt, id, "folder")
  #
  # dt_z <- cbind(dt_z,dt)

  # browser()

  #
  # meta_learner_fit <- meta_learner$fit(dt_z)

  # learners on the full dataset
  # full_train_list <- lapply(learners, function(f) f$fit(dt))
  #
  # step_0_predictions <- mapply(
  #   function(f, model, newdata)
  #     f$predictor(model = model, newdata = newdata),
  #   learners,
  #   full_train_list,
  #   MoreArgs = list(newdata = dt)
  # )
  #
  # step_0_predictions<- apply(as.matrix(step_0_predictions),
  #                  MARGIN = 2,
  #                  log)

  # Name the columns
  # colnames(step_0_predictions) <- z_covariates
  #
  # setDT(cbind(as.data.frame.matrix(step_0_predictions),
  #             dt[,c("node",
  #                   "tij")]
  #             ))
  #
  # fitted_values <- meta_learner$predictor(meta_learner_fit,
  #                                   dt_z)


  out <- list(
    learners = learners,
    metalearner = meta_learner,
    superlearner=meta_learner_fits,
    data_info = list(
      id = id,
      status = status,
      event_time = event_time,
      nodes = nodes,
      nfold = nfold
    )
  )

  class(out) <- "poisson_superlearner"

  return(out)



}






