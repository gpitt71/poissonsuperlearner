#' Individual Data Pre-Processing
#'
#' This function predicts the survivor function or the hazard using the superlearner.
#'
#' @return predictions
#'
#' @export
predict.poisson_superlearner <- function(object,
                                         newdata,
                                         type = "survival",
                                         k=1,
                                         ...) {


  setDT(newdata)
  dt <- copy(newdata)

  z_covariates <- paste0("Z", 1:length(object$learners))


  # Predict on the validation set your pseudo-observations ----
  learners_predictions <- mapply(
    function(f, model, newdata)
      f$predictor(model = model, newdata = dt),
    object$learners,
    object$superlearner[[k]]$learners_fit,
    MoreArgs = list(newdata = dt)
  )

  pseudo_observations_data <- apply(as.matrix(learners_predictions), MARGIN = 2, log)

  learners_sf <- apply(as.matrix(learners_predictions), MARGIN = 2, function(x) exp(-cumsum(x)) )

  colnames(learners_sf) <- paste0("survival_function_l",1:ncol(learners_sf))

  # Name the columns

  colnames(pseudo_observations_data) <- z_covariates

  setDT(as.data.frame.matrix(pseudo_observations_data))

  dt_pred <- object$superlearner[[k]]$model$predictor(object$superlearner[[k]]$meta_learner_fit, newdata =
                                                  cbind(pseudo_observations_data,dt))




  dt[['pwch_times_tij']] <- dt_pred
  dt[['survival_function']] <- exp(-cumsum(dt_pred))


  dt <- cbind(dt,learners_sf)


  return(dt)




}
