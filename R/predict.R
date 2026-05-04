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

  n_crisks <- object$data_info$n_crisks

  if (length(cause) != 1L || is.na(cause) || cause != as.integer(cause)) {
    stop("'cause' must be a single positive integer.", call. = FALSE)
  }

  cause <- as.integer(cause)

  if (cause < 1L || cause > n_crisks) {
    stop(
      sprintf("'cause' must be between 1 and %d.", n_crisks),
      call. = FALSE
    )
  }

  ## ------------------------------------------------------------
  ## Cause-specific object components
  ## ------------------------------------------------------------

  learners_by_cause <- object$learners_by_cause
  learners_labels_by_cause <- object$learners_labels_by_cause
  z_covariates_by_cause <- object$z_covariates_by_cause
  metalearner_by_cause <- object$metalearner_by_cause

  ## Backward compatibility for older fitted objects
  if (is.null(learners_by_cause)) {
    if (is.null(object$learners)) {
      stop("No learner library found in the fitted object.", call. = FALSE)
    }

    learners_by_cause <- replicate(
      n_crisks,
      object$learners,
      simplify = FALSE
    )
  }

  if (is.null(learners_labels_by_cause)) {
    learners_labels_by_cause <- object$data_info$learners_labels_by_cause
  }

  if (is.null(learners_labels_by_cause)) {
    learners_labels_by_cause <- lapply(learners_by_cause, function(x) {
      labs <- names(x)
      if (is.null(labs)) {
        labs <- paste0("learner_", seq_along(x))
      }
      labs
    })
  }

  if (is.null(z_covariates_by_cause)) {
    z_covariates_by_cause <- object$data_info$z_covariates_by_cause
  }

  if (is.null(z_covariates_by_cause)) {
    z_covariates_by_cause <- lapply(
      learners_by_cause,
      function(x) paste0("Z", seq_along(x))
    )
  }

  if (is.null(metalearner_by_cause)) {
    metalearner_by_cause <- replicate(
      n_crisks,
      object$metalearner,
      simplify = FALSE
    )
  }

  n_learners_by_cause <- lengths(learners_by_cause)

  if (length(learners_by_cause) != n_crisks) {
    stop("'object$learners_by_cause' must have length equal to n_crisks.", call. = FALSE)
  }

  if (length(learners_labels_by_cause) != n_crisks) {
    stop("'object$learners_labels_by_cause' must have length equal to n_crisks.", call. = FALSE)
  }

  if (length(z_covariates_by_cause) != n_crisks) {
    stop("'object$z_covariates_by_cause' must have length equal to n_crisks.", call. = FALSE)
  }

  if (any(n_learners_by_cause == 0L)) {
    stop("At least one cause has no retained learners.", call. = FALSE)
  }


  ## ------------------------------------------------------------
  ## Resolve model selector cause-by-cause
  ## ------------------------------------------------------------

  resolve_model_by_cause <- function(model) {

    if (is.null(model) || length(model) != 1L) {
      stop("'model' must be a scalar.", call. = FALSE)
    }

    if (is.numeric(model)) {
      if (is.na(model) || model != as.integer(model)) {
        stop("Numeric 'model' must be one of 0, 1, 2, ...", call. = FALSE)
      }

      model <- as.integer(model)

      if (model == 0L) {
        return(list(
          type = "sl",
          index_by_cause = rep(0L, n_crisks),
          label = "sl"
        ))
      }

      bad <- which(model > n_learners_by_cause)

      if (length(bad) > 0L) {
        stop(
          sprintf(
            "Numeric model %d is unavailable for cause(s): %s.",
            model,
            paste(bad, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      return(list(
        type = "learner",
        index_by_cause = rep(model, n_crisks),
        label = paste0("learner_", model)
      ))
    }

    if (!is.character(model) || is.na(model)) {
      stop("'model' must be a character scalar or a numeric scalar.", call. = FALSE)
    }

    model_chr <- trimws(model)
    model_chr_lc <- tolower(model_chr)

    if (model_chr_lc %in% c("sl", "superlearner", "super_learner")) {
      return(list(
        type = "sl",
        index_by_cause = rep(0L, n_crisks),
        label = "sl"
      ))
    }

    if (grepl("^learner_[0-9]+$", model_chr_lc)) {
      j <- as.integer(sub("^learner_", "", model_chr_lc))
      bad <- which(j > n_learners_by_cause)

      if (length(bad) > 0L) {
        stop(
          sprintf(
            "'%s' is unavailable for cause(s): %s.",
            model_chr,
            paste(bad, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      return(list(
        type = "learner",
        index_by_cause = rep(j, n_crisks),
        label = model_chr
      ))
    }

    index_by_cause <- integer(n_crisks)

    for (k in seq_len(n_crisks)) {
      labels_k <- learners_labels_by_cause[[k]]
      j <- match(model_chr, labels_k)

      if (is.na(j)) {
        stop(
          sprintf(
            "Learner label '%s' is not available for cause %d.",
            model_chr,
            k
          ),
          call. = FALSE
        )
      }

      index_by_cause[k] <- j
    }

    list(
      type = "learner",
      index_by_cause = index_by_cause,
      label = model_chr
    )
  }

  model_sel <- resolve_model_by_cause(model)


  ## ------------------------------------------------------------
  ## Build prediction data
  ## ------------------------------------------------------------

  data.table::setDT(newdata)

  tmp <- data.table::copy(newdata)
  tmp[, internal_psl_ix := seq_len(.N)]

  tmp <- tmp[
    ,
    setdiff(
      names(tmp),
      c(
        object$data_info$event_time,
        object$data_info$status,
        object$data_info$id
      )
    ),
    with = FALSE
  ]

  cond_zero <- 0 %in% times
  all_zero <- all(times == 0)
  cond_times_larger_than_max <- times > object$data_info$maximum_followup

  pwch_cols <- paste0("pwch_", seq_len(n_crisks))

  if (all(cond_times_larger_than_max)) {
    warning(
      paste0(
        "All the entries in the input times are larger than the maximum follow-up: ",
        as.character(object$data_info$maximum_followup)
      )
    )
    return(NULL)
  }

  vec_dt <- data.table::data.table(
    tmp_prediction_time = times[times <= object$data_info$maximum_followup]
  )

  data.table::setnames(
    vec_dt,
    "tmp_prediction_time",
    object$data_info$event_time
  )

  tmp[, dummy := 1L]
  vec_dt[, dummy := 1L]

  data_pp <- merge(
    tmp,
    vec_dt,
    by = "dummy",
    allow.cartesian = TRUE
  )[, dummy := NULL]

  if (is.null(data_pp[[object$data_info$id]])) {
    data_pp[[object$data_info$id]] <- seq_len(nrow(data_pp))
  }

  if (is.null(data_pp[[object$data_info$status]])) {
    data_pp[[object$data_info$status]] <- 0
  }


  ## ------------------------------------------------------------
  ## Handle time zero
  ## ------------------------------------------------------------

  zero_time <- NULL

  if (cond_zero) {
    tmptcol <- object$data_info$event_time

    zero_time <- data_pp[get(tmptcol) == 0]

    for (cl in pwch_cols) {
      data.table::set(zero_time, j = cl, value = 0)
    }

    data.table::set(zero_time, j = "survival_function", value = 1)
    data.table::set(zero_time, j = "absolute_risk", value = 0)

    data_pp <- data_pp[get(tmptcol) != 0]
  }

  if (!all_zero) {

    data_pp <- data_pre_processing(
      data_pp,
      id = object$data_info$id,
      status = object$data_info$status,
      predictions = TRUE,
      event_time = object$data_info$event_time,
      nodes = object$data_info$nodes
    )

    if (!is.null(object$data_info$variable_transformation)) {
      apply_transformations(
        data_pp,
        object$data_info$variable_transformation
      )
    }

    data_pp[, deltatime := tij]
    data_pp[, tij := 1]


    ## ------------------------------------------------------------
    ## Cause-specific prediction
    ## ------------------------------------------------------------

    dt_pred <- vector("list", n_crisks)

    for (k in seq_len(n_crisks)) {

      learners_k <- learners_by_cause[[k]]
      z_k <- z_covariates_by_cause[[k]]
      sl_k <- object$superlearner[[k]]
      n_k <- n_learners_by_cause[k]

      use_sl_k <- isTRUE(model_sel$type == "sl") &&
        n_k > 1L &&
        !is.null(sl_k$meta_learner_fit) &&
        !is.null(metalearner_by_cause[[k]])

      if (use_sl_k) {

        learners_predictions_k <- lapply(seq_len(n_k), function(j) {
          learners_k[[j]]$private_predictor(
            model = sl_k$learners_fit[[j]],
            newdata = data_pp
          )
        })

        learners_predictions_k <- do.call(
          cbind,
          lapply(learners_predictions_k, as.numeric)
        )

        learners_predictions_k <- log(learners_predictions_k)

        colnames(learners_predictions_k) <- z_k

        pseudo_observations_k <- data.table::as.data.table(
          learners_predictions_k
        )

        dt_pred[[k]] <- metalearner_by_cause[[k]]$private_predictor(
          model = sl_k$meta_learner_fit,
          newdata = cbind(pseudo_observations_k, data_pp)
        )

      } else {

        learner_index_k <- if (model_sel$type == "learner") {
          model_sel$index_by_cause[k]
        } else {
          1L
        }

        fit_k <- if (n_k == 1L) {
          sl_k$learners_fit
        } else {
          sl_k$learners_fit[[learner_index_k]]
        }

        dt_pred[[k]] <- learners_k[[learner_index_k]]$private_predictor(
          model = fit_k,
          newdata = data_pp
        )
      }
    }

    data_pp[, (pwch_cols) := dt_pred]

    sum_of_hazards <- paste(pwch_cols, collapse = " + ")
    eval(parse(
      text = paste0("data_pp[, pwch_dot := ", sum_of_hazards, "]")
    ))

    mapply(
      function(pwch, name) {
        data_pp[
          ,
          (paste0("cumulative_hazard_", name)) :=
            cumsum(get(pwch) * deltatime),
          by = id
        ]
      },
      pwch_cols,
      gsub("pwch_", "", pwch_cols)
    )

    haz <- as.matrix(
      data_pp[, .SD, .SDcols = patterns("^pwch_[0-9]+$")]
    )

    S <- pch_survival(
      id = data_pp$id,
      dt = data_pp$deltatime,
      haz = haz
    )

    data_pp[, survival_function := S]

    data_pp[
      ,
      absolute_risk := pch_absolute_risk(
        id,
        deltatime,
        haz,
        cause_idx = cause
      )
    ]

    data_pp <- data_pp[
      ,
      .SD[.N],
      by = id
    ][
      ,
      times := node_start + deltatime
    ]

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

  } else {

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

    d <- zero_time[, ..columns_ss]
  }


  ## ------------------------------------------------------------
  ## Add time zero rows, if needed
  ## ------------------------------------------------------------

  if (cond_zero && !all_zero) {
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

    zero_time <- zero_time[, ..columns_ss]
    d <- rbind(zero_time, d, fill = TRUE)
  }


  ## ------------------------------------------------------------
  ## Add rows for times larger than maximum follow-up
  ## ------------------------------------------------------------

  if (any(cond_times_larger_than_max)) {

    vec_dt2 <- data.table::data.table(
      tmp_prediction_time = times[cond_times_larger_than_max]
    )

    data.table::setnames(
      vec_dt2,
      "tmp_prediction_time",
      object$data_info$event_time
    )

    tmp2 <- data.table::copy(tmp)
    tmp2[, dummy := 1L]
    vec_dt2[, dummy := 1L]

    d2 <- merge(
      tmp2,
      vec_dt2,
      by = "dummy",
      allow.cartesian = TRUE
    )[, dummy := NULL]

    d2[, c(pwch_cols, "survival_function", "absolute_risk") := NA_real_]

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

    d2 <- d2[, ..columns_ss]

    d <- rbind(d, d2, fill = TRUE)
  }


  ## ------------------------------------------------------------
  ## Restore user-facing id
  ## ------------------------------------------------------------

  if (object$data_info$id %in% names(d)) {
    d[, (object$data_info$id) := NULL]
  }

  data.table::setnames(
    d,
    old = "internal_psl_ix",
    new = object$data_info$id
  )

  d <- d[order(get(object$data_info$id), get(object$data_info$event_time))]

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
  # data_pp <- data_pp[, .SD[.N], by = id][, times := as.numeric(as.character(node)) +
  #                                          deltatime]
  data_pp <- data_pp[, .SD[.N], by = id][, times := node_start + deltatime]


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


  if (object$data_info$id %in% names(d)) {
    d[, (object$data_info$id) := NULL]
  }
  setnames(d, new = object$data_info$id, old = "internal_psl_ix")
  d <- d[order(get(object$data_info$id))]
  return(d)




}
