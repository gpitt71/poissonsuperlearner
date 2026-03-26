#' Predict hazards, survival and absolute risk from a fitted Poisson Super Learner
#'
#' Computes **cause-specific piecewise-constant hazards** (`pwch_k`), the corresponding
#' **survival function**, and **absolute risk** for a given cause, at user-supplied
#' prediction horizons `times`, for each row in `newdata`.
#'
#' Internally, `newdata` is expanded to a Cartesian product with the requested
#' `times`, converted to long Poisson format on `object$data_info$nodes`, and hazards
#' are predicted either from the stacked super learner (`model = "sl"`) or from one
#' selected fitted base learner. Survival and absolute risk are then computed from
#' the predicted hazards.
#'
#' @param object `poisson_superlearner`. A fitted ensemble from [Superlearner()].
#' @param newdata `data.frame`/`data.table`. New covariate data (one row per subject).
#'   If `newdata` contains the original `event_time`, `status`, or `id` columns used
#'   for fitting, they are ignored for prediction.
#' @param times `numeric`. Prediction horizon(s). May include `0`.
#'   Times larger than `object$data_info$maximum_followup` are not supported:
#'   if **all** requested times exceed the maximum follow-up, a warning is issued and
#'   `NULL` is returned; if only **some** exceed, output rows for those times are
#'   returned with `NA` predictions.
#' @param cause `numeric(1)`. Cause index (1, 2, ...) used for the `absolute_risk`
#'   calculation.
#' @param model Scalar model selector. Default is `"sl"` for the stacked super learner.
#'   Other allowed values are:
#'   \describe{
#'     \item{`0` or `"sl"`}{Use the super learner prediction.}
#'     \item{learner label}{Use one stored base learner by its label in
#'       `object$data_info$learners_labels`.}
#'     \item{`"learner_j"`}{Use the `j`-th stored learner.}
#'     \item{integer `j >= 1`}{Use the `j`-th stored learner.}
#'   }
#'   Numeric positions refer to the learners actually stored in the fitted object.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' **Special case `times = 0`:** when `0` is included in `times`, the returned rows
#' have `survival_function = 1`, `absolute_risk = 0`, and all `pwch_k = 0` at time 0.
#'
#' **Identifiers in the output:** if `newdata` contains the `id` column, it is carried
#' into the output. If `newdata` does not contain an id column, an internal id is
#' created for computation, but it is not guaranteed to appear in the returned table
#' unless it was present in `newdata`.
#'
#' @return A `data.table` with one row per `(row in newdata, time in times)` and columns:
#' \describe{
#'   \item{(original columns)}{All columns from `newdata` (excluding ignored event columns).}
#'   \item{time column}{A column with name `object$data_info$event_time` holding the requested horizon.}
#'   \item{pwch_1, pwch_2, ...}{Predicted cause-specific piecewise hazards at the horizon.}
#'   \item{survival_function}{Predicted survival probability at the horizon.}
#'   \item{absolute_risk}{Predicted cumulative incidence (absolute risk) for `cause` at the horizon.}
#' }
#'
#' @examples
#' d <- simulateStenoT1(30, competing_risks = TRUE)
#'
#' learners <- list(
#'   lasso = Learner_glmnet(
#'     covariates = "sex",
#'     alpha = 1,
#'     lambda = 0.01,
#'     cross_validation = FALSE
#'   ),
#'   ridge = Learner_glmnet(
#'     covariates = c("sex", "value_LDL"),
#'     alpha = 0,
#'     lambda = 0.01,
#'     cross_validation = FALSE
#'   )
#' )
#'
#' fit <- Superlearner(
#'   data = d,
#'   id = "id",
#'   status = "status_cvd",
#'   event_time = "time_cvd",
#'   learners = learners,
#'   number_of_nodes = 3,
#'   nfold = 2
#' )
#' p <- predict(fit, newdata = d[1:3], times = c(0, 2), cause = 1)
#' p[, .(id, time_cvd, absolute_risk)]
#'
#'
#' @export
predict.poisson_superlearner <- function(object,
                                         newdata,
                                         times,
                                         cause = 1,
                                         model = "sl",
                                         ...) {

  model_sel <- resolve_prediction_model(object, model)

  setDT(newdata)
  tmp <- copy(newdata)
  tmp[, internal_psl_ix := 1:.N]

  tmp <- tmp[, setdiff(names(tmp), c(object$data_info$event_time,
                                     object$data_info$status,
                                     object$data_info$id)), with = FALSE]

  cond_zero <- 0 %in% times
  all_zero <- all(times == 0)
  cond_times_larger_than_max <- times > object$data_info$maximum_followup

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
  } else {
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
    data_pp <- merge(tmp, vec_dt, by = "dummy", allow.cartesian = TRUE)[, dummy := NULL]
  }

  if (is.null(data_pp[[object$data_info$id]])) {
    data_pp[[object$data_info$id]] <- 1:nrow(data_pp)
  }

  if (is.null(data_pp[[object$data_info$status]])) {
    data_pp[[object$data_info$status]] <- 0
  }

  if (cond_zero) {
    tmptcol <- object$data_info$event_time
    zero_time <- data_pp[get(tmptcol) == 0]
    for (cl in pwch_cols) set(zero_time, j = cl, value = 0)
    set(zero_time, j = "survival_function", value = 1)
    set(zero_time, j = "absolute_risk", value = 0)
    data_pp <- data_pp[get(tmptcol) != 0]
  }

  if (all_zero) {
    return(zero_time)
  }

  data_pp <- data_pre_processing(
    data_pp,
    id = object$data_info$id,
    status = object$data_info$status,
    predictions = TRUE,
    event_time = object$data_info$event_time,
    nodes = object$data_info$nodes
  )

  if (!is.null(object$data_info$variable_transformation)) {
    apply_transformations(data_pp, object$data_info$variable_transformation)
  }

  z_covariates <- paste0("Z", seq_along(object$learners))

  data_pp[, deltatime := tij][, tij := 1]

  use_superlearner <- isTRUE(model_sel$type == "sl") &&
    length(object$learners) > 1L &&
    !is.null(object$metalearner)

  if (use_superlearner) {

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
      colnames(out) <- z_covariates
      as.data.table(as.data.frame.matrix(out))
    })

    dt_pred <- mapply(
      function(superlearner, pseudo_observations_data) {
        object$metalearner$private_predictor(
          model = superlearner$meta_learner_fit,
          newdata = cbind(pseudo_observations_data, data_pp)
        )
      },
      object$superlearner,
      pseudo_observations_data,
      SIMPLIFY = FALSE
    )

  } else {

    learner_index <- if (model_sel$type == "learner") {
      model_sel$index
    } else {
      1L
    }

    dt_pred <- lapply(seq_len(object$data_info$n_crisks), function(k) {
      fit_k <- if (length(object$learners) == 1L) {
        object$superlearner[[k]]$learners_fit
      } else {
        object$superlearner[[k]]$learners_fit[[learner_index]]
      }

      object$learners[[learner_index]]$private_predictor(
        model = fit_k,
        newdata = data_pp
      )
    })
  }

  data_pp[, paste0("pwch_", 1:object$data_info$n_crisks) := dt_pred]

  sum_of_hazards <- paste(pwch_cols, collapse = " + ")
  pwch_dot_string <- paste0("data_pp[, pwch_dot :=", sum_of_hazards, "]")
  eval(parse(text = pwch_dot_string))

  mapply(function(pwch, name) {
    data_pp[, (paste0("cumulative_hazard_", name)) := cumsum(get(pwch) * deltatime), by = id]
  }, pwch_cols, gsub("pwch_", "", pwch_cols))

  haz <- as.matrix(data_pp[, .SD, .SDcols = patterns("^pwch_[0-9]+$")])
  S <- pch_survival(id = data_pp$id, dt = data_pp$deltatime, haz = haz)
  data_pp[, survival_function := S]

  data_pp[, absolute_risk := pch_absolute_risk(id, deltatime, haz, cause_idx = cause)]

  data_pp <- data_pp[, .SD[.N], by = id][, times := as.numeric(as.character(node)) + deltatime]

  columns_ss <- unique(
    c(
      colnames(newdata),
      object$data_info$event_time,
      pwch_cols,
      "survival_function",
      "absolute_risk",
      "internal_psl_ix"
    )
  )

  d <- data_pp[, ..columns_ss]

  if (cond_zero) {
    d <- rbind(zero_time, d)
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
    d2[, c(pwch_cols, "survival_function", "absolute_risk") := NA_real_]

    if (object$data_info$id %in% colnames(d)) {
      d2[[object$data_info$id]] <- (nrow(data_pp) + 1):(nrow(data_pp) + nrow(d2))
    }

    d <- rbind(d, d2, fill = TRUE)
  }

  d[, (object$data_info$id) := NULL]
  setnames(d, new = object$data_info$id, old = "internal_psl_ix")
  d <- d[order(get(object$data_info$id))]
  return(d)
}



#' Predict hazards, survival and absolute risk from a fitted base learner
#'
#' Computes **cause-specific piecewise-constant hazards** (`pwch_k`), the corresponding
#' **survival function**, and **absolute risk** for a given cause, at user-supplied
#' prediction horizons `times`, using a fitted `base_learner` object (single learner;
#' no stacking).
#'
#' Internally, `newdata` is expanded to a Cartesian product with `times`, converted to
#' long Poisson format on `object$data_info$nodes`, and the fitted learner for each
#' cause in `object$learner_fit` is used to predict the cause-specific hazards.
#' Survival and absolute risk are then computed from the predicted hazards.
#'
#' @param object `base_learner`. A fitted object returned by [fit_learner()].
#'   It contains the learner specification in `object$model` and cause-specific fitted
#'   models in `object$learner_fit`.
#' @param newdata `data.frame`/`data.table`. New covariate data (one row per subject).
#'   If `newdata` contains the original `event_time`, `status`, or `id` columns used
#'   for fitting, they are ignored for prediction.
#' @param times `numeric`. Prediction horizon(s). May include `0`.
#'   Times larger than `object$data_info$maximum_followup` are not supported:
#'   if **all** requested times exceed the maximum follow-up, a warning is issued and
#'   `NULL` is returned; if only **some** exceed, output rows for those times are
#'   returned with `NA` predictions.
#' @param cause `numeric(1)`. Cause index (1, 2, ...) used for the `absolute_risk`
#'   calculation.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' **Special case `times = 0`:** when `0` is included in `times`, the returned rows
#' have `survival_function = 1`, `absolute_risk = 0`, and all `pwch_k = 0` at time 0.
#'
#' **Identifiers in the output:** if `newdata` contains the `id` column, it is carried
#' into the output. If `newdata` does not contain an id column, an internal id is
#' created for computation, but it is not guaranteed to appear in the returned table
#' unless it was present in `newdata`.
#'
#' @return A `data.table` with one row per `(row in newdata, time in times)` and columns:
#' \describe{
#'   \item{(original columns)}{All columns from `newdata` (excluding ignored event columns).}
#'   \item{time column}{A column with name `object$data_info$event_time` holding the requested horizon.}
#'   \item{pwch_1, pwch_2, ...}{Predicted cause-specific piecewise hazards at the horizon.}
#'   \item{survival_function}{Predicted survival probability at the horizon.}
#'   \item{absolute_risk}{Predicted cumulative incidence (absolute risk) for `cause` at the horizon.}
#' }
#'
#' @examples
#' d <- simulateStenoT1(120, competing_risks = TRUE)
#' lrn <- Learner_glmnet(covariates = c("age", "value_LDL"), lambda = 0, cross_validation = FALSE)
#' bl <- fit_learner(d, learner = lrn, id="id", status="status_cvd", event_time="time_cvd",
#'                   number_of_nodes=8)
#' p <- predict(bl, newdata = d[1:5], times = c(0, 2, 5), cause = 1)
#' head(p)
#'
#' @export
predict.base_learner <- function(object,
                                         newdata,
                                         times,
                                         cause = 1,
                                         ...) {

  setDT(newdata)

  tmp <- copy(newdata)
  tmp[,internal_psl_ix:=1:.N]
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




# no problem writing over id
if (is.null(data_pp[[object$data_info$id]])) {
  data_pp[[object$data_info$id]] <- 1:nrow(data_pp)
}

if (is.null(data_pp[[object$data_info$status]])) {
  data_pp[[object$data_info$status]] <- 0
}


  if(cond_zero){


    tmptcol <- object$data_info$event_time  # character vector
    zero_time <- data_pp[get(tmptcol) == 0]
    for (cl in pwch_cols) set(zero_time, j = cl, value = 0)
    set(zero_time, j = "survival_function", value = 1)
    set(zero_time, j = "absolute_risk", value = 0)
    data_pp<-data_pp[get(tmptcol) != 0]
  }

  if(all_zero){

    return(zero_time)
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
      "absolute_risk",
      "internal_psl_ix"
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
    d2[, c(pwch_cols, "survival_function", "absolute_risk") := NA_real_]


    if (object$data_info$id %in% colnames(d)) {
      d2[[object$data_info$id]] <- (nrow(data_pp) + 1):(nrow(data_pp) + nrow(d2))
    }



    d <- rbind(d, d2, fill = TRUE)


  }


  d[, (object$data_info$id) := NULL]
  setnames(d, new = object$data_info$id, old = "internal_psl_ix")
  d <- d[order(get(object$data_info$id))]
  return(d)




}
