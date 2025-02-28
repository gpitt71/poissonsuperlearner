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
                         status="status",
                         event_time="time",
                         learners,
                         nodes=NULL,
                         nfold = 3){



  # save some relevant values
  n <- nrow(data)
  id_fold <- sample(1:nfold,n,replace=TRUE,prob = rep(1/nfold,nfold))

  dt_id <- data.table(
    folder = id_fold,
    id = data[[id]]
  )

  # pre-process the data
  dt <- data_pre_processing(
    data = data,
    id = id,
    status = status,
    nodes = nodes,
    event_time = event_time
  )

  dt <- merge(dt,
              dt_id,
              by="id",
              all.x = T)



  dt_z <- NULL
  z_covariates <- paste0("Z", 1:length(learners))
  # Train your models in the training set ----
  for(ix in 1:nfold){

    # browser()

    tmp_train <- dt[folder != ix,]
    tmp_val <- dt[folder == ix,]

    train_list <- lapply(learners, function(f) f$fit(tmp_train))

  # Predict on the validation set your pseudo-observations ----
    val_list <- mapply(
      function(f, model, newdata)
        f$predictor(model = model, newdata = newdata),
      learners,
      train_list,
      MoreArgs = list(newdata = tmp_val)
    )

    val_list<- apply(as.matrix(val_list),
    MARGIN = 2,
    log)

    # Name the columns

    colnames(val_list) <- z_covariates

    dt_z <- rbind(dt_z,
                  data.table(val_list)[,c("id",
                                          "folder"):=list(tmp_val$id,
                                                          ix)])


  }


  # browser()

  setorder(dt_z, id, "folder")
  setorder(dt, id, "folder")

  dt_z <- cbind(dt_z,dt)


  # We do another round of glmnet for combining the predictors ----
  combiner <- Learner_glm(covariates = z_covariates,
                          cross_validation=FALSE)


  superlearner <- combiner$fit(dt_z)

  out <- list(
    learners = learners,
    superlearner=superlearner,
    nodes=nodes,
    nfold=nfold)

  class(out) <- "poisson_superlearner"

  return(out)



}






