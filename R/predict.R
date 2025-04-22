#' Individual Data Pre-Processing
#'
#' This function predicts the survivor function or the hazard using the superlearner.
#'
#' @return predictions
#'
#' @export
predict.poisson_superlearner <- function(object, newdata, type = "survival", k=1,step=.5, ...) {
  # We think about this one.. Mostly dataformatting issues
  # pre-process the data
  # dt <- data_pre_processing(
  #   data = newdata,
  #   id = object$data_info[['id']],
  #   status = object$data_info[['status']],
  #   nodes = object$data_info[['nodes']],
  #   event_time = object$data_info[['event_time']]
  # )
  # browser()
  dt <- newdata
  setDT(dt)

  # train_list <- lapply(object$superlearner$learners_fit, function(f) f$fit(dt))
  #train_list <- lapply(object$learners, function(f) f$fit(tmp_05))
  z_covariates <- paste0("Z", 1:length(object$learners))


  # Predict on the validation set your pseudo-observations ----
  learners_predictions <- mapply(
    function(f, model, newdata)
      f$predictor(model = model, newdata = newdata),
    object$learners,
    object$superlearner[[k]]$learners_fit,
    MoreArgs = list(newdata = dt)
  )

  pseudo_observations_data <- apply(as.matrix(learners_predictions), MARGIN = 2, log)

  # Name the columns

  colnames(pseudo_observations_data) <- z_covariates

  setDT(as.data.frame.matrix(pseudo_observations_data))


  dt_pred <- object$superlearner[[1]]$model$predictor(object$superlearner[[1]]$meta_learner_fit, newdata =
                                                  cbind(pseudo_observations_data,dt))




  dt[['pw_constant_hazard']] <- dt_pred

  dt[, survival_function := exp(-cumsum(pw_constant_hazard * step))]

  return(dt)




}
