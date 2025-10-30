#' Compute Poisson deviance for Poisson Super Learner objects
#'
#' This method evaluates the out-of-sample Poisson deviance of the fitted
#' learners that form the ensemble as well as the meta learner on new data.
#'
#' @param object A fitted object of class `poisson_superlearner` returned by
#'   [Superlearner()].
#' @param ... Additional arguments passed to methods. Currently unused.
#'
#' @return A list with two elements: `learners`, a `data.table` containing the
#'   Poisson deviance for each base learner, and `metalearner`, a `data.table`
#'   with the deviance for the meta learner. If the super learner was trained
#'   with a single learner, the `metalearner` entry is `NULL`.
#'
#' @export
poisson_deviance <- function(object, ...) {
  UseMethod("poisson_deviance")
}

#' @rdname poisson_deviance
#' @param newdata A `data.frame` containing the data that should be used to
#'   evaluate the Poisson deviance.
#' @export
poisson_deviance.poisson_superlearner <- function(object, newdata, ...) {
  browser()
  if (missing(newdata)) {
    stop("Argument 'newdata' must be supplied.")
  }

  if (!inherits(object, "poisson_superlearner")) {
    stop("Input object must be of class 'poisson_superlearner'.")
  }

  newdata <- copy(as.data.table(newdata))

  id_col <- object$data_info$id
  status_col <- object$data_info$status
  event_time_col <- object$data_info$event_time
  start_time_col <- object$data_info$start_time
  end_time_col <- object$data_info$end_time

  if (is.null(newdata[[id_col]])) {
    newdata[, (id_col) := seq_len(.N)]
  }

  if (is.null(newdata[[status_col]])) {
    stop("The status column specified during training is missing in 'newdata'.")
  }

  if (isTRUE(object$data_info$interval_data_type)) {
    if (is.null(newdata[[start_time_col]]) || is.null(newdata[[end_time_col]])) {
      stop("The start and end time columns are required for interval data.")
    }

    dt <- data_pre_processing_interval_data(
      data = newdata,
      id = id_col,
      status = status_col,
      start_time = start_time_col,
      end_time = end_time_col,
      nodes = object$data_info$nodes
    )
  } else {
    if (is.null(event_time_col) || is.null(newdata[[event_time_col]])) {
      stop("The event time column specified during training is missing in 'newdata'.")
    }

    uncensored_01 <- !(0 %in% newdata[[status_col]])

    dt <- data_pre_processing(
      data = newdata,
      id = id_col,
      status = status_col,
      event_time = event_time_col,
      nodes = object$data_info$nodes,
      uncensored_01 = uncensored_01
    )
  }

  if (!is.null(object$data_info$variable_transformation)) {
    apply_transformations(dt, object$data_info$variable_transformation)
  }



  data_by_competing_risk <- split(dt, by = "k")


  if (length(object$learners) == 1L) {



    dt_pred <- lapply(object$superlearner, function(sl_fit) {
      object$learners[[1]]$private_predictor(
        model = sl_fit$learners_fit, newdata = dt[k==1,])
    })




  } else{
    #

    learners_predictions <- mapply(
      function(crisk_cause,
               superlearner,
               newdata,
               learners)
        learners_hat(crisk_cause, superlearner, newdata, learners),
      crisk_cause = as.list(1:object$data_info$n_crisks),
      superlearner = object$superlearner,
      MoreArgs = list(newdata = data_pp, learners = object$learners),
      SIMPLIFY = FALSE
    )


    pseudo_observations_data <- lapply(learners_predictions, function(mx) {
      out <- matrix(
        apply(as.matrix(
          mx,
          nrow = nrow(newdata),
          ncol = length(z_covariates)
        ), MARGIN = 2, log),
        nrow = nrow(data_pp),
        ncol = length(z_covariates)
      )
      # out <-as.matrix(mx)
      colnames(out) <- z_covariates # Name the columns
      as.data.table(as.data.frame.matrix(out))
    })


    dt_pred <- mapply(
      function(crisk_cause,
               superlearner,
               pseudo_observations_data) {
        superlearner$model$predictor(superlearner$meta_learner_fit,
                                     newdata = cbind(pseudo_observations_data, data_pp))
      },
      as.list(1:object$data_info$n_crisks),
      object$superlearner,
      pseudo_observations_data,
      SIMPLIFY = F
    )


  }



  if (!"folder" %in% names(dt)) {
    dt[, folder := 1L]
  }

  learner_names <- names(object$learners)
  if (is.null(learner_names) || any(!nzchar(learner_names))) {
    learner_names <- paste0("learner_", seq_along(object$learners))
  }

  n_crisks <- object$data_info$n_crisks
  z_covariates <- paste0("Z", seq_along(object$learners))

  # Identify the name of the selected meta learner (if any)
  metalearner_name <- NULL
  if (!is.null(object$metalearner)) {
    if (!is.null(object$meta_learner_cross_validation)) {
      candidate_meta <- setdiff(
        unique(object$meta_learner_cross_validation$model),
        learner_names
      )
      if (length(candidate_meta) == 1L) {
        metalearner_name <- candidate_meta
      } else if (length(candidate_meta) > 1L) {
        metalearner_name <- object$meta_learner_cross_validation[
          model %in% candidate_meta
        ][which.min(deviance), model]
      }
    }
    if (is.null(metalearner_name) || length(metalearner_name) == 0L) {
      metalearner_name <- "metalearner"
    }
  }

  # Helper to compute Poisson deviance contributions
  dev_contribution <- function(delta, mu) {
    out <- numeric(length(delta))
    positive <- delta > 0

    out[!positive] <- mu[!positive]

    valid_mu <- mu > 0 & positive
    out[valid_mu] <- delta[valid_mu] * log(delta[valid_mu] / mu[valid_mu]) -
      (delta[valid_mu] - mu[valid_mu])

    # observations with events but zero/negative mean -> infinite deviance
    out[positive & !valid_mu] <- Inf

    # for zero counts with zero/negative mean, contribution is zero
    out[!positive & mu <= 0] <- 0

    out
  }

  # Compute predictions for learners and meta learner
  learner_long <- vector("list", n_crisks)
  meta_long <- if (!is.null(object$metalearner)) vector("list", n_crisks) else NULL

  for (cause in seq_len(n_crisks)) {
    dt_cause <- dt[k == cause]
    if (nrow(dt_cause) == 0L) {
      learner_long[[cause]] <- NULL
      if (!is.null(meta_long)) meta_long[[cause]] <- NULL
      next
    }

    preds_matrix <- learners_hat(
      crisk_cause = cause,
      superlearner = object$superlearner[[cause]],
      newdata = dt_cause,
      learners = object$learners
    )

    preds_matrix <- as.matrix(preds_matrix)
    if (ncol(preds_matrix) == 1L && length(learner_names) == 1L) {
      preds_matrix <- matrix(preds_matrix, ncol = 1L)
    }

    preds_dt <- as.data.table(preds_matrix)
    setnames(preds_dt, learner_names)

    preds_dt[, `:=`(
      id = dt_cause$id,
      node = dt_cause$node,
      folder = dt_cause$folder,
      delta = dt_cause$deltaij,
      tij = dt_cause$tij,
      cause = cause
    )]

    learner_long[[cause]] <- melt(
      preds_dt,
      id.vars = c("id", "node", "folder", "delta", "tij", "cause"),
      variable.name = "model",
      value.name = "lambda",
      variable.factor = FALSE
    )

    if (!is.null(object$metalearner)) {
      pseudo_mat <- log(as.matrix(preds_dt[, ..learner_names]))
      colnames(pseudo_mat) <- z_covariates

      meta_input <- cbind(as.data.table(pseudo_mat), dt_cause)
      meta_lambda <- object$superlearner[[cause]]$model$predictor(
        object$superlearner[[cause]]$meta_learner_fit,
        newdata = meta_input
      )

      meta_long[[cause]] <- data.table(
        model = metalearner_name,
        id = dt_cause$id,
        folder = dt_cause$folder,
        delta = dt_cause$deltaij,
        tij = dt_cause$tij,
        lambda = meta_lambda
      )
    }
  }

  learner_long <- rbindlist(learner_long, use.names = TRUE, fill = TRUE)

  compute_deviance <- function(dt_long) {
    if (nrow(dt_long) == 0L) {
      return(data.table(model = character(), deviance = numeric()))
    }

    dt_long <- dt_long[is.finite(lambda) & !is.na(tij) & !is.na(delta)]
    if (nrow(dt_long) == 0L) {
      return(data.table(model = character(), deviance = numeric()))
    }

    mu <- dt_long$lambda * dt_long$tij
    dt_long[, contribution := dev_contribution(delta, mu)]

    dev_id <- dt_long[, .(deviance_i = sum(contribution)), by = .(model, folder, id)]
    dev_folder <- dev_id[, .(deviance_v = 2 * sum(deviance_i)), by = .(model, folder)]
    dev_folder[, .(deviance = mean(deviance_v)), by = model]
  }

  learners_deviance <- compute_deviance(learner_long)

  metalearner_deviance <- NULL
  if (!is.null(object$metalearner)) {
    meta_long <- rbindlist(meta_long, use.names = TRUE, fill = TRUE)
    metalearner_deviance <- compute_deviance(meta_long)
  }

  list(
    learners = learners_deviance,
    metalearner = metalearner_deviance
  )
}
