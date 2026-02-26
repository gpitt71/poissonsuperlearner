#' Poisson Super-Learner Predictions
#'
#' This method computes the survival function and the absolute risk prediction of a \code{poisson_superlearner} object, based on a given data set, at given \code{times} for a given \code{cause}.
#'
#'
#' @param object \code{poisson_superlearner} for absolute risk and survival function predictions.
#' @param newdata \code{data.frame}, new data to predict the absolute risk and the survival function for.
#' @param times \code{numeric}, time(s) at which to predict the absolute risk and survival function.
#' @param cause \code{numeric}, competing risk to predict the absolute risk and the survival function for.
#'
#' @return \code{data.table} containing for each row of \code{newdata} the \code{poisson_superlearner} predictions for the survival function and the absolute risk predictions for some \code{cause} at the given \code{times}.
#'
#' @export
predict.poisson_superlearner <- function(object,
                                         newdata,
                                         times,
                                         cause = 1,
                                         ...) {

  setDT(newdata)

  tmp <- copy(newdata)
  # here we disregard the event_time column if present in the newdata
  tmp <- tmp[, setdiff(names(tmp), c(object$data_info$event_time,
                                     object$data_info$status,
                                     object$data_info$id)), with = FALSE]


  ## checks on the data
  cond_zero <- 0 %in% times

  all_zero <- all(times==0)

  cond_times_larger_than_max <- times > object$data_info$maximum_followup

  ## frame hazard problem
  pwch_cols <- paste0("pwch_", 1:object$data_info$n_crisks)

  if (all(cond_times_larger_than_max)) {
    warning(
      paste0(
        "All the entries in the input times are larger than the maximum follow-up: ",
        as.character(object$data_info$maximum_followup)
      )
    )
    d <- NULL

    return(d)

  } else{
    eval(parse(
      text = paste0(
        "
    vec_dt <- data.table(

    ",
        object$data_info$event_time,
        " = times[times <= object$data_info$maximum_followup]
  )
    "
      )
    ))

    tmp[, dummy := 1]
    vec_dt[, dummy := 1]

    # Merge on dummy to create Cartesian product
    data_pp <- merge(tmp, vec_dt, by = "dummy", allow.cartesian = TRUE)[, dummy := NULL]
  }
  # }



  if(cond_zero){

    zero_time <- data_pp[time==0,]
    zero_time[, c(pwch_cols, 'survival_function', 'absolute_risk') := list(rep(0, length(pwch_cols)), 1,0)]
    data_pp<-data_pp[time!=0,]
  }

  if(all_zero){

    return(zero_time)
  }

  # no problem writing over id
  if (is.null(data_pp[[object$data_info$id]])) {
    data_pp[[object$data_info$id]] <- 1:nrow(data_pp)
  }

  if (is.null(data_pp[[object$data_info$status]])) {
    data_pp[[object$data_info$status]] <- 0
  }

  data_pp <- data_pre_processing(
    data_pp,
    id = object$data_info$id,
    status = object$data_info$status,
    predictions=TRUE,
    event_time = object$data_info$event_time,
    nodes = object$data_info$nodes
  )


  # In case of variable transformations ----
  if (!is.null(object$data_info$variable_transformation)) {apply_transformations(data_pp, object$data_info$variable_transformation)}


  # Set covariates for metalearner
  z_covariates <- paste0("Z", 1:length(object$learners))


  # Predict on the validation set your pseudo-observations ----
  #

  data_pp[, deltatime := tij][, tij := 1]

  if (length(object$learners) == 1L) {

    dt_pred <- lapply(object$superlearner, function(sl_fit) {
      object$learners[[1]]$private_predictor(
        model = sl_fit$learners_fit, newdata = data_pp)
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
        superlearner$model$private_predictor(superlearner$meta_learner_fit,
                                     newdata = cbind(pseudo_observations_data, data_pp))
      },
      as.list(1:object$data_info$n_crisks),
      object$superlearner,
      pseudo_observations_data,
      SIMPLIFY = F
    )


  }

  #

  # save casue-specific pwch

  data_pp[, paste0("pwch_", 1:object$data_info$n_crisks) := dt_pred]



  # save sum of pwch

  sum_of_hazards <- paste(pwch_cols, collapse = " + ")

  pwch_dot_string <- paste0("data_pp[, pwch_dot :=", sum_of_hazards, "]")

  eval(parse(text = pwch_dot_string))


  # compute cumulative hazard

  mapply(function(pwch, name) {
    data_pp[, (paste0("cumulative_hazard_", name)) := cumsum(get(pwch) * deltatime), by = id]
  }, pwch_cols, gsub("pwch_", "", pwch_cols))


  # compute survival function

  ## c++
  haz <- as.matrix(data_pp[, .SD, .SDcols = patterns("^pwch_[0-9]+$")])
  S <- pch_survival(id = data_pp$id, dt = data_pp$deltatime, haz = haz)
  data_pp[, survival_function := S]


  data_pp[, absolute_risk := pch_absolute_risk(id, deltatime, haz, cause_idx = cause)]

  ####
  data_pp <- data_pp[, .SD[.N], by = id][, times := as.numeric(as.character(node)) +
                                           deltatime]



  columns_ss <- unique(
    c(
      colnames(newdata),
      object$data_info$event_time,
      pwch_cols,
      "survival_function",
      "absolute_risk"
    )
  )

  d <- data_pp[, ..columns_ss]


  if (cond_zero) {

    d<- rbind(zero_time,
              d)

  }



  if (any(cond_times_larger_than_max)) {
    eval(parse(
      text = paste0(
        "
    vec_dt2 <- data.table(

    ",
        object$data_info$event_time,
        " = times[cond_times_larger_than_max]
  )
    "
      )
    ))

    tmp[, dummy := 1]
    vec_dt2[, dummy := 1]

    d2 <- merge(tmp, vec_dt2, by = "dummy", allow.cartesian = TRUE)[, dummy := NULL]
    d2[, c(pwch_cols, 'survival_function') := list(rep(NA, length(pwc_cols)), NA)]


    if (object$data_info$id %in% colnames(d)) {
      d2[[object$data_info$id]] <- (nrow(data_pp) + 1):(nrow(data_pp) + nrow(d2))
    }



    d <- rbind(d, d2)


  }


  return(d)




}




#' Learner Predictions
#'
#' This method computes the survival function and the absolute risk prediction of a \code{base_learner} object, based on a given data set, at given \code{times} for a given \code{cause}.
#'
#'
#' @param object \code{poisson_superlearner} for absolute risk and survival function predictions.
#' @param newdata \code{data.frame}, new data to predict the absolute risk and the survival function for.
#' @param times \code{numeric}, time(s) at which to predict the absolute risk and survival function.
#' @param cause \code{numeric}, competing risk to predict the absolute risk and the survival function for.
#'
#' @return \code{data.table} containing for each row of \code{newdata} the \code{poisson_superlearner} predictions for the survival function and the absolute risk predictions for some \code{cause} at the given \code{times}.
#'
#' @export
predict.base_learner <- function(object,
                                         newdata,
                                         times,
                                         cause = 1,
                                         ...) {

  setDT(newdata)

  tmp <- copy(newdata)
  # here we disregard the event_time column if present in the newdata
  tmp <- tmp[, setdiff(names(tmp), c(object$data_info$event_time,
                                     object$data_info$status,
                                     object$data_info$id)), with = FALSE]


  ## checks on the data
  cond_zero <- 0 %in% times

  all_zero <- all(times==0)

  cond_times_larger_than_max <- times > object$data_info$maximum_followup

  ## frame hazard problem
  pwch_cols <- paste0("pwch_", 1:object$data_info$n_crisks)

  if (all(cond_times_larger_than_max)) {
    warning(
      paste0(
        "All the entries in the input times are larger than the maximum follow-up: ",
        as.character(object$data_info$maximum_followup)
      )
    )
    d <- NULL

    return(d)

  } else{
    eval(parse(
      text = paste0(
        "
    vec_dt <- data.table(

    ",
        object$data_info$event_time,
        " = times[times <= object$data_info$maximum_followup]
  )
    "
      )
    ))

    # # no problem writing over id
    # if (is.null(tmp[[object$data_info$id]])) {
    #   tmp[[object$data_info$id]] <- 1:nrow(tmp)
    # }
    #
    # if (is.null(tmp[[object$data_info$status]])) {
    #   tmp[[object$data_info$status]] <- 0
    # }

    tmp[, dummy := 1]
    vec_dt[, dummy := 1]

    # Merge on dummy to create Cartesian product
    data_pp <- merge(tmp, vec_dt, by = "dummy", allow.cartesian = TRUE)[, dummy := NULL]
  }
  # }



  if(cond_zero){

    zero_time <- data_pp[time==0,]
    zero_time[, c(pwch_cols, 'survival_function', 'absolute_risk') := list(rep(0, length(pwch_cols)), 1,0)]
    data_pp<-data_pp[time!=0,]
  }

  if(all_zero){

    return(zero_time)
  }

  # no problem writing over id
  if (is.null(data_pp[[object$data_info$id]])) {
    data_pp[[object$data_info$id]] <- 1:nrow(data_pp)
  }

  if (is.null(data_pp[[object$data_info$status]])) {
    data_pp[[object$data_info$status]] <- 0
  }

  data_pp <- data_pre_processing(
    data_pp,
    id = object$data_info$id,
    status = object$data_info$status,
    predictions=TRUE,
    event_time = object$data_info$event_time,
    nodes = object$data_info$nodes
  )



  #  ()
  if (!is.null(object$data_info$variable_transformation)) {

    apply_transformations(data_pp, object$data_info$variable_transformation)
  }


  # Set covariates for metalearner
  z_covariates <- paste0("Z", 1:length(object$learners))


  # Predict on the validation set your pseudo-observations ----
  #

  data_pp[, deltatime := tij][, tij := 1]


    dt_pred <- lapply(object$learner_fit, function(x) {
      object$model$private_predictor(
        model = x, newdata = data_pp)
    })


  # save casue-specific pwch

  data_pp[, paste0("pwch_", 1:object$data_info$n_crisks) := dt_pred]



  # save sum of pwch

  sum_of_hazards <- paste(pwch_cols, collapse = " + ")

  pwch_dot_string <- paste0("data_pp[, pwch_dot :=", sum_of_hazards, "]")

  eval(parse(text = pwch_dot_string))


  # compute cumulative hazard

  mapply(function(pwch, name) {
    data_pp[, (paste0("cumulative_hazard_", name)) := cumsum(get(pwch) * deltatime), by = id]
  }, pwch_cols, gsub("pwch_", "", pwch_cols))


  # compute survival function

  ## c++
  haz <- as.matrix(data_pp[, .SD, .SDcols = patterns("^pwch_[0-9]+$")])
  S <- pch_survival(id = data_pp$id, dt = data_pp$deltatime, haz = haz)
  data_pp[, survival_function := S]


  data_pp[, absolute_risk := pch_absolute_risk(id, deltatime, haz, cause_idx = cause)]

  # abs_risk approx
  # data_pp[, survival_function_shift := shift(survival_function, fill = 1), by =
  #           id]
  # absolute_risk_string <- paste0(
  #   "data_pp[, absolute_risk_2 := cumsum(survival_function_shift * pwch_",
  #   cause,
  #   "*deltatime), by = id]"
  # )
  # eval(parse(text = absolute_risk_string))


  ## non c++
  # hazard_terms <- paste0("cumulative_hazard_", 1:object$data_info$n_crisks)
  # sum_expr <- paste(pwch_cols, collapse = " + ")
  # survival_function_string <- paste0("data_pp[, survival_function := exp(-cumsum((", sum_expr, ")*deltatime)),by=id]")
  # eval(parse(text = survival_function_string))

  # shift survival function
  # data_pp[, survival_function_shift := shift(survival_function, fill = 1), by =
  #           id]
  # absolute_risk_string <- paste0(
  #   "data_pp[, absolute_risk := cumsum(survival_function_shift * pwch_",
  #   cause,
  #   "/pwch_dot * (1-exp(-pwch_dot*deltatime))), by = id]"
  # )
  # eval(parse(text = absolute_risk_string))

  ####
  data_pp <- data_pp[, .SD[.N], by = id][, times := as.numeric(as.character(node)) +
                                           deltatime]



  columns_ss <- unique(
    c(
      colnames(newdata),
      object$data_info$event_time,
      pwch_cols,
      "survival_function",
      "absolute_risk"
    )
  )

  d <- data_pp[, ..columns_ss]


  if (cond_zero) {

    d<- rbind(zero_time,
              d)

  }



  if (any(cond_times_larger_than_max)) {
    eval(parse(
      text = paste0(
        "
    vec_dt2 <- data.table(

    ",
        object$data_info$event_time,
        " = times[cond_times_larger_than_max]
  )
    "
      )
    ))

    tmp[, dummy := 1]
    vec_dt2[, dummy := 1]

    d2 <- merge(tmp, vec_dt2, by = "dummy", allow.cartesian = TRUE)[, dummy := NULL]
    d2[, c(pwch_cols, 'survival_function') := list(rep(NA, length(pwc_cols)), NA)]


    if (object$data_info$id %in% colnames(d)) {
      d2[[object$data_info$id]] <- (nrow(data_pp) + 1):(nrow(data_pp) + nrow(d2))
    }



    d <- rbind(d, d2)


  }


  return(d)




}
