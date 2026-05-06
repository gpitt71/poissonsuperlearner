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

  id_col <- object$data_info$id
  status_col <- object$data_info$status
  event_time_col <- object$data_info$event_time
  grid_nodes <- object$data_info$nodes
  maximum_followup <- object$data_info$maximum_followup
  n_crisks <- object$data_info$n_crisks

  times <- as.numeric(times)

  if (!length(times)) {
    stop("`times` must contain at least one value.", call. = FALSE)
  }

  if (anyNA(times)) {
    stop("`times` cannot contain NA values.", call. = FALSE)
  }

  if (any(times < 0)) {
    stop("`times` must be non-negative.", call. = FALSE)
  }

  if (length(cause) != 1L || cause < 1L || cause > n_crisks) {
    stop(
      "`cause` must be a single integer between 1 and object$data_info$n_crisks.",
      call. = FALSE
    )
  }

  learners_by_cause <- object$learners

  if (
    is.null(learners_by_cause) ||
    !is.list(learners_by_cause) ||
    length(learners_by_cause) != n_crisks
  ) {
    if (!is.null(object$learners_by_cause)) {
      learners_by_cause <- object$learners_by_cause
    }
  }

  if (
    is.null(learners_by_cause) ||
    !is.list(learners_by_cause) ||
    length(learners_by_cause) != n_crisks
  ) {
    stop(
      "Could not find cause-specific learners in `object$learners`.",
      call. = FALSE
    )
  }

  pwch_cols <- paste0("pwch_", seq_len(n_crisks))

  learner_labels <- function(k) {

    labs <- NULL

    if (
      !is.null(object$data_info$learners_labels) &&
      is.list(object$data_info$learners_labels) &&
      length(object$data_info$learners_labels) >= k
    ) {
      labs <- object$data_info$learners_labels[[k]]
    }

    if (
      is.null(labs) ||
      length(labs) != length(learners_by_cause[[k]]) ||
      anyNA(labs) ||
      any(!nzchar(labs))
    ) {
      labs <- names(learners_by_cause[[k]])
    }

    if (
      is.null(labs) ||
      length(labs) != length(learners_by_cause[[k]]) ||
      anyNA(labs) ||
      any(!nzchar(labs))
    ) {
      labs <- paste0("learner_", seq_along(learners_by_cause[[k]]))
    }

    labs
  }
  best_discrete_sl_index_by_cause <- function() {

    cv_dt <- object$cross_validation_deviance

    index_by_cause <- integer(n_crisks)
    label_by_cause <- character(n_crisks)

    for (k in seq_len(n_crisks)) {

      labs_k <- learner_labels(k)
      n_k <- length(labs_k)

      ## If only one learner is retained for this cause, it is automatically discrete SL.
      if (n_k == 1L) {
        index_by_cause[k] <- 1L
        label_by_cause[k] <- labs_k[1L]
        next
      }

      if (
        is.null(cv_dt) ||
        !is.data.frame(cv_dt) ||
        nrow(cv_dt) == 0L
      ) {
        stop(
          "`model = 'discrete_sl'` requires `object$cross_validation_deviance`, unless each cause has only one retained learner.",
          call. = FALSE
        )
      }

      cv_df <- as.data.frame(cv_dt)

      if (!("deviance" %in% names(cv_df))) {
        stop(
          "`object$cross_validation_deviance` must contain a `deviance` column for `model = 'discrete_sl'`.",
          call. = FALSE
        )
      }

      if ("cause_index" %in% names(cv_df)) {
        cv_k <- cv_df[as.integer(cv_df[["cause_index"]]) == k, , drop = FALSE]
      } else if (
        "cause" %in% names(cv_df) &&
        !is.null(names(learners_by_cause)) &&
        length(names(learners_by_cause)) >= k &&
        nzchar(names(learners_by_cause)[k])
      ) {
        cv_k <- cv_df[as.character(cv_df[["cause"]]) == names(learners_by_cause)[k], , drop = FALSE]
      } else {
        stop(
          "`object$cross_validation_deviance` must contain `cause_index` or cause labels matching `object$learners` for `model = 'discrete_sl'`.",
          call. = FALSE
        )
      }

      if (nrow(cv_k) == 0L) {
        stop(
          sprintf(
            "No cross-validation deviance is available for cause %d, so `model = 'discrete_sl'` cannot choose among %d learners.",
            k,
            n_k
          ),
          call. = FALSE
        )
      }

      dev_k <- as.numeric(cv_k[["deviance"]])

      if (all(is.na(dev_k) | !is.finite(dev_k))) {
        stop(
          sprintf(
            "All cross-validation deviances are missing or non-finite for cause %d.",
            k
          ),
          call. = FALSE
        )
      }

      dev_k[is.na(dev_k) | !is.finite(dev_k)] <- Inf
      best_row <- which.min(dev_k)

      ix <- NA_integer_

      ## Prefer label matching, because it is robust to pruning/re-indexing.
      if ("learner" %in% names(cv_k)) {
        best_label <- as.character(cv_k[["learner"]][best_row])
        ix <- match(best_label, labs_k)
      }

      ## Fall back to learner_index if labels are unavailable or do not match.
      if (
        is.na(ix) &&
        "learner_index" %in% names(cv_k)
      ) {
        ix <- as.integer(cv_k[["learner_index"]][best_row])
      }

      if (is.na(ix) || ix < 1L || ix > n_k) {
        stop(
          sprintf(
            "The best cross-validation learner for cause %d could not be matched to the retained learner library. Available learners: %s.",
            k,
            paste(sprintf("%d='%s'", seq_along(labs_k), labs_k), collapse = ", ")
          ),
          call. = FALSE
        )
      }

      index_by_cause[k] <- ix
      label_by_cause[k] <- labs_k[ix]
    }

    list(
      type = "learner",
      index_by_cause = index_by_cause,
      label_by_cause = label_by_cause
    )
  }
  resolve_model <- function(model) {

    if (is.null(model) || !length(model)) {
      stop(
        "`model` must be 'sl', a learner label, a learner index, or a vector of cause-specific learner labels/indices.",
        call. = FALSE
      )
    }

    if (!(length(model) %in% c(1L, n_crisks))) {
      stop(
        sprintf(
          "`model` must have length 1 or length %d, one entry per cause.",
          n_crisks
        ),
        call. = FALSE
      )
    }

    ## ------------------------------------------------------------
    ## Numeric selector
    ## ------------------------------------------------------------

    if (is.numeric(model) || is.integer(model)) {

      if (anyNA(model) || any(model != as.integer(model))) {
        stop("Numeric `model` values must be integer learner indices.", call. = FALSE)
      }

      model_int <- as.integer(model)

      if (length(model_int) == 1L && model_int == 0L) {
        return(list(
          type = "sl",
          index_by_cause = rep(0L, n_crisks),
          label_by_cause = rep("sl", n_crisks)
        ))
      }

      if (any(model_int < 1L)) {
        stop(
          "Numeric `model` values must be positive learner indices. Use scalar 0 or 'sl' for the Super Learner.",
          call. = FALSE
        )
      }

      index_by_cause <- if (length(model_int) == 1L) {
        rep(model_int, n_crisks)
      } else {
        model_int
      }

      for (k in seq_len(n_crisks)) {

        labs_k <- learner_labels(k)

        if (index_by_cause[k] > length(labs_k)) {
          stop(
            sprintf(
              "Learner index %d is not available for cause %d. Available learners: %s.",
              index_by_cause[k],
              k,
              paste(sprintf("%d='%s'", seq_along(labs_k), labs_k), collapse = ", ")
            ),
            call. = FALSE
          )
        }
      }

      return(list(
        type = "learner",
        index_by_cause = index_by_cause,
        label_by_cause = vapply(
          seq_len(n_crisks),
          function(k) learner_labels(k)[index_by_cause[k]],
          character(1L)
        )
      ))
    }

    ## ------------------------------------------------------------
    ## Character selector
    ## ------------------------------------------------------------

    if (!is.character(model)) {
      stop(
        "`model` must be 'sl', a learner label, a learner index, or a vector of cause-specific learner labels/indices.",
        call. = FALSE
      )
    }

    if (anyNA(model) || any(!nzchar(trimws(model)))) {
      stop("Character `model` values cannot be missing or empty.", call. = FALSE)
    }

    model_chr <- trimws(model)
    model_lc <- tolower(model_chr)

    sl_values <- c("sl", "superlearner", "super_learner")

    discrete_sl_values <- c(
      "discrete_sl",
      "discrete_superlearner",
      "discrete_super_learner",
      "discretesl"
    )

    if (
      length(model_chr) == 1L &&
      model_lc %in% sl_values
    ) {
      return(list(
        type = "sl",
        index_by_cause = rep(0L, n_crisks),
        label_by_cause = rep("sl", n_crisks)
      ))
    }

    if (
      length(model_chr) == 1L &&
      model_lc %in% discrete_sl_values
    ) {
      return(best_discrete_sl_index_by_cause())
    }

    if (any(model_lc %in% c(sl_values, discrete_sl_values))) {
      stop(
        "`model = 'sl'` and `model = 'discrete_sl'` must be supplied as scalar model selectors. Do not mix them with cause-specific base-learner selectors.",
        call. = FALSE
      )
    }

    model_chr_by_cause <- if (length(model_chr) == 1L) {
      rep(model_chr, n_crisks)
    } else {
      model_chr
    }

    index_by_cause <- integer(n_crisks)
    label_by_cause <- character(n_crisks)

    for (k in seq_len(n_crisks)) {

      labs_k <- learner_labels(k)
      selector_k <- model_chr_by_cause[k]
      selector_k_lc <- tolower(selector_k)

      if (grepl("^[0-9]+$", selector_k)) {

        ix <- as.integer(selector_k)

      } else if (grepl("^learner_[0-9]+$", selector_k_lc)) {

        ix <- as.integer(sub("^learner_", "", selector_k_lc))

      } else {

        ix <- match(selector_k, labs_k)
      }

      if (is.na(ix) || ix < 1L || ix > length(labs_k)) {
        stop(
          sprintf(
            "Learner '%s' is not available for cause %d. Available learners: %s.",
            selector_k,
            k,
            paste(sprintf("%d='%s'", seq_along(labs_k), labs_k), collapse = ", ")
          ),
          call. = FALSE
        )
      }

      index_by_cause[k] <- ix
      label_by_cause[k] <- labs_k[ix]
    }

    list(
      type = "learner",
      index_by_cause = index_by_cause,
      label_by_cause = label_by_cause
    )
  }

  model_sel <- resolve_model(model)

  learner_index_for_cause <- function(k) {
    model_sel$index_by_cause[k]
  }

  get_fit <- function(k, m) {

    fits_k <- object$superlearner[[k]]$learners_fit

    if (length(learners_by_cause[[k]]) == 1L) {
      return(fits_k)
    }

    fits_k[[m]]
  }

  predict_base <- function(k, m, endpoint_data, n_expanded) {

    pred <- learners_by_cause[[k]][[m]]$private_predictor(
      model = get_fit(k, m),
      newdata = endpoint_data,
      grid_nodes = grid_nodes
    )

    if (length(pred) != n_expanded) {
      stop(
        sprintf(
          "Prediction length mismatch for cause %s, learner %s: got %s, expected %s.",
          k,
          m,
          length(pred),
          n_expanded
        ),
        call. = FALSE
      )
    }

    as.numeric(pred)
  }

  predict_meta_glmnet_hazard <- function(meta_fit, z_mat, deltatime) {

    eps <- 1e-15

    z_log <- log(pmax(z_mat, eps))
    storage.mode(z_log) <- "double"

    deltatime <- as.numeric(deltatime)

    out <- rep(NA_real_, nrow(z_log))

    ok <- stats::complete.cases(z_log) &
      is.finite(deltatime) &
      deltatime > 0

    if (any(ok)) {

      mu <- stats::predict(
        meta_fit,
        newx = z_log[ok, , drop = FALSE],
        newoffset = log(deltatime[ok]),
        s = 0,
        type = "response"
      )

      ## glmnet returns the Poisson mean:
      ##   E[N_ij] = deltatime_ij * hazard_ij
      ## The rest of the prediction code needs the hazard.
      out[ok] <- as.numeric(mu) / deltatime[ok]
    }

    out
  }

  build_endpoint_data <- function(combo_dt) {

    grid_dt <- data.table::data.table(node = grid_nodes)
    grid_dt[, prev_node_psl := data.table::shift(node)]
    grid_dt[, (event_time_col) := node]

    endpoint_dt <- grid_dt[
      combo_dt,
      on = event_time_col,
      roll = Inf
    ]

    endpoint_dt[
      ,
      node_start := data.table::fifelse(
        get(event_time_col) == node,
        prev_node_psl,
        node
      )
    ]

    endpoint_dt <- endpoint_dt[!is.na(node_start)]

    endpoint_dt[
      ,
      c("node", "tij") := list(
        node_start,
        get(event_time_col) - node_start
      )
    ]

    endpoint_dt[, deltaij := 0L]

    endpoint_dt[]
  }


  ## ------------------------------------------------------------
  ## One row for each newdata x times combination
  ## ------------------------------------------------------------

  newdata_dt <- data.table::as.data.table(data.table::copy(newdata))
  newdata_dt[, .psl_internal_ix := seq_len(.N)]

  tmp <- data.table::copy(newdata_dt)

  drop_cols <- intersect(
    c(id_col, status_col, event_time_col),
    names(tmp)
  )

  if (length(drop_cols)) {
    tmp[, (drop_cols) := NULL]
  }

  n_new <- nrow(tmp)
  n_times <- length(times)

  combo_dt <- tmp[rep(seq_len(n_new), each = n_times)]
  combo_dt[, .psl_time_ix := rep(seq_along(times), times = n_new)]
  combo_dt[, .psl_pred_id := seq_len(.N)]
  combo_dt[, (event_time_col) := rep(times, times = n_new)]

  zero_combo <- combo_dt[get(event_time_col) == 0]
  future_combo <- combo_dt[get(event_time_col) > maximum_followup]
  valid_combo <- combo_dt[
    get(event_time_col) > 0 &
      get(event_time_col) <= maximum_followup
  ]

  if (nrow(valid_combo) == 0L && nrow(zero_combo) == 0L) {
    warning(
      paste0(
        "All entries in `times` are larger than the maximum follow-up: ",
        as.character(maximum_followup),
        "."
      ),
      call. = FALSE
    )
    return(NULL)
  }


  ## ------------------------------------------------------------
  ## t = 0
  ## ------------------------------------------------------------

  zero_out <- NULL

  if (nrow(zero_combo)) {
    zero_out <- data.table::copy(zero_combo)

    for (cc in pwch_cols) {
      zero_out[, (cc) := 0]
    }

    zero_out[, survival_function := 1]
    zero_out[, absolute_risk := 0]
  }


  ## ------------------------------------------------------------
  ## t > maximum follow-up
  ## ------------------------------------------------------------

  future_out <- NULL

  if (nrow(future_combo)) {
    future_out <- data.table::copy(future_combo)

    for (cc in pwch_cols) {
      future_out[, (cc) := NA_real_]
    }

    future_out[, survival_function := NA_real_]
    future_out[, absolute_risk := NA_real_]
  }


  ## ------------------------------------------------------------
  ## 0 < t <= maximum follow-up
  ## ------------------------------------------------------------

  valid_out <- NULL

  if (nrow(valid_combo)) {

    endpoint_data <- build_endpoint_data(valid_combo)

    if (!is.null(object$data_info$variable_transformation)) {
      apply_transformations(
        endpoint_data,
        object$data_info$variable_transformation
      )
    }

    ## This is the important part:
    ## expand each requested subject-time endpoint into the Poisson
    ## representation over all intervals up to that time.
    data_pp <- make_validation_skeleton(
      valid_data = endpoint_data,
      grid_nodes = grid_nodes,
      cause = 1L,
      node_col = "node",
      tij_col = "tij",
      event_col = "deltaij",
      keep_cols = ".psl_pred_id"
    )

    if (!(".psl_pred_id" %in% names(data_pp))) {
      stop(
        "`make_validation_skeleton()` must preserve `.psl_pred_id` during prediction.",
        call. = FALSE
      )
    }

    if (!("tij" %in% names(data_pp))) {
      stop(
        "`make_validation_skeleton()` must return `tij`.",
        call. = FALSE
      )
    }

    data_pp[, deltatime := as.numeric(tij)]
    data_pp <- data_pp[deltatime > 0]

    if ("node" %in% names(data_pp)) {
      data.table::setorderv(data_pp, c(".psl_pred_id", "node"))
    } else {
      data.table::setorderv(data_pp, ".psl_pred_id")
    }

    n_expanded <- nrow(data_pp)

    dt_pred <- vector("list", n_crisks)

    for (k in seq_len(n_crisks)) {

      n_learners_k <- length(learners_by_cause[[k]])
      meta_fit_k <- object$superlearner[[k]]$meta_learner_fit

      use_meta_k <- isTRUE(model_sel$type == "sl") &&
        n_learners_k > 1L &&
        !is.null(meta_fit_k)

      if (use_meta_k) {

        z_mat <- matrix(
          NA_real_,
          nrow = n_expanded,
          ncol = n_learners_k
        )

        colnames(z_mat) <- paste0("Z", seq_len(n_learners_k))

        for (m in seq_len(n_learners_k)) {
          z_mat[, m] <- predict_base(
            k = k,
            m = m,
            endpoint_data = endpoint_data,
            n_expanded = n_expanded
          )
        }

        dt_pred[[k]] <- predict_meta_glmnet_hazard(
          meta_fit = meta_fit_k,
          z_mat = z_mat,
          deltatime = data_pp[["deltatime"]]
        )

      } else {

        learner_ix <- if (isTRUE(model_sel$type == "learner")) {
          learner_index_for_cause(k)
        } else {
          1L
        }

        dt_pred[[k]] <- predict_base(
          k = k,
          m = learner_ix,
          endpoint_data = endpoint_data,
          n_expanded = n_expanded
        )
      }
    }

    data_pp[, (pwch_cols) := dt_pred]

    haz <- as.matrix(data_pp[, .SD, .SDcols = pwch_cols])

    data_pp[
      ,
      survival_function := pch_survival(
        id = .psl_pred_id,
        dt = deltatime,
        haz = haz
      )
    ]

    data_pp[
      ,
      absolute_risk := pch_absolute_risk(
        .psl_pred_id,
        deltatime,
        haz,
        cause_idx = cause
      )
    ]

    last_pp <- data_pp[
      ,
      .SD[.N],
      by = .psl_pred_id
    ]

    pred_cols <- c(
      ".psl_pred_id",
      pwch_cols,
      "survival_function",
      "absolute_risk"
    )

    valid_out <- merge(
      valid_combo,
      last_pp[, ..pred_cols],
      by = ".psl_pred_id",
      all.x = TRUE,
      sort = FALSE
    )
  }


  ## ------------------------------------------------------------
  ## Assemble output
  ## ------------------------------------------------------------

  out <- data.table::rbindlist(
    list(zero_out, valid_out, future_out),
    use.names = TRUE,
    fill = TRUE
  )

  data.table::setorderv(
    out,
    c(".psl_internal_ix", ".psl_time_ix")
  )

  data.table::setnames(
    out,
    old = ".psl_internal_ix",
    new = id_col
  )

  private_cols <- c(".psl_pred_id", ".psl_time_ix")

  out[, (intersect(private_cols, names(out))) := NULL]

  return_cols <- unique(
    c(
      id_col,
      setdiff(
        names(newdata_dt),
        c(id_col, status_col, event_time_col, ".psl_internal_ix")
      ),
      event_time_col,
      pwch_cols,
      "survival_function",
      "absolute_risk"
    )
  )

  return_cols <- intersect(return_cols, names(out))

  out <- out[, ..return_cols]

  out[]
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

  id_col <- object$data_info$id
  status_col <- object$data_info$status
  event_time_col <- object$data_info$event_time
  grid_nodes <- object$data_info$nodes
  maximum_followup <- object$data_info$maximum_followup
  n_crisks <- object$data_info$n_crisks

  times <- as.numeric(times)

  if (!length(times)) {
    stop("`times` must contain at least one value.", call. = FALSE)
  }

  if (anyNA(times)) {
    stop("`times` cannot contain NA values.", call. = FALSE)
  }

  if (any(times < 0)) {
    stop("`times` must be non-negative.", call. = FALSE)
  }

  if (
    length(cause) != 1L ||
    is.na(cause) ||
    cause < 1L ||
    cause > n_crisks
  ) {
    stop(
      "`cause` must be a single integer between 1 and object$data_info$n_crisks.",
      call. = FALSE
    )
  }

  if (
    is.null(object$learner_fit) ||
    !is.list(object$learner_fit) ||
    length(object$learner_fit) != n_crisks
  ) {
    stop(
      "`object$learner_fit` must be a list with one fitted model per cause.",
      call. = FALSE
    )
  }

  if (
    is.null(object$model) ||
    !is.function(object$model$private_predictor)
  ) {
    stop(
      "`object$model` must contain a learner with private_predictor().",
      call. = FALSE
    )
  }

  pwch_cols <- paste0("pwch_", seq_len(n_crisks))


  ## ------------------------------------------------------------
  ## One row for each newdata x times combination
  ## ------------------------------------------------------------

  newdata_dt <- data.table::as.data.table(data.table::copy(newdata))
  newdata_dt[, .psl_internal_ix := seq_len(.N)]

  tmp <- data.table::copy(newdata_dt)

  drop_cols <- intersect(
    c(id_col, status_col, event_time_col),
    names(tmp)
  )

  if (length(drop_cols)) {
    tmp[, (drop_cols) := NULL]
  }

  n_new <- nrow(tmp)
  n_times <- length(times)

  combo_dt <- tmp[rep(seq_len(n_new), each = n_times)]
  combo_dt[, .psl_time_ix := rep(seq_along(times), times = n_new)]
  combo_dt[, .psl_pred_id := seq_len(.N)]
  combo_dt[, (event_time_col) := rep(times, times = n_new)]

  zero_combo <- combo_dt[get(event_time_col) == 0]
  future_combo <- combo_dt[get(event_time_col) > maximum_followup]
  valid_combo <- combo_dt[
    get(event_time_col) > 0 &
      get(event_time_col) <= maximum_followup
  ]

  if (nrow(valid_combo) == 0L && nrow(zero_combo) == 0L) {
    warning(
      paste0(
        "All entries in `times` are larger than the maximum follow-up: ",
        as.character(maximum_followup),
        "."
      ),
      call. = FALSE
    )

    return(NULL)
  }


  ## ------------------------------------------------------------
  ## t = 0
  ## ------------------------------------------------------------

  zero_out <- NULL

  if (nrow(zero_combo)) {
    zero_out <- data.table::copy(zero_combo)

    for (cc in pwch_cols) {
      zero_out[, (cc) := 0]
    }

    zero_out[, survival_function := 1]
    zero_out[, absolute_risk := 0]
  }


  ## ------------------------------------------------------------
  ## t > maximum follow-up
  ## ------------------------------------------------------------

  future_out <- NULL

  if (nrow(future_combo)) {
    future_out <- data.table::copy(future_combo)

    for (cc in pwch_cols) {
      future_out[, (cc) := NA_real_]
    }

    future_out[, survival_function := NA_real_]
    future_out[, absolute_risk := NA_real_]
  }


  ## ------------------------------------------------------------
  ## 0 < t <= maximum follow-up
  ## ------------------------------------------------------------

  valid_out <- NULL

  if (nrow(valid_combo)) {

    grid_dt <- data.table::data.table(node = grid_nodes)
    grid_dt[, prev_node_psl := data.table::shift(node)]
    grid_dt[, (event_time_col) := node]

    endpoint_data <- grid_dt[
      valid_combo,
      on = event_time_col,
      roll = Inf
    ]

    endpoint_data[
      ,
      node_start := data.table::fifelse(
        get(event_time_col) == node,
        prev_node_psl,
        node
      )
    ]

    endpoint_data <- endpoint_data[!is.na(node_start)]

    endpoint_data[
      ,
      c("node", "tij") := list(
        node_start,
        get(event_time_col) - node_start
      )
    ]

    endpoint_data[, deltaij := 0L]

    if (!is.null(object$data_info$variable_transformation)) {
      apply_transformations(
        endpoint_data,
        object$data_info$variable_transformation
      )
    }

    data_pp <- make_validation_skeleton(
      valid_data = endpoint_data,
      grid_nodes = grid_nodes,
      cause = 1L,
      node_col = "node",
      tij_col = "tij",
      event_col = "deltaij",
      keep_cols = ".psl_pred_id"
    )

    if (!(".psl_pred_id" %in% names(data_pp))) {
      stop(
        "`make_validation_skeleton()` must preserve `.psl_pred_id` during prediction.",
        call. = FALSE
      )
    }

    if (!("tij" %in% names(data_pp))) {
      stop(
        "`make_validation_skeleton()` must return `tij`.",
        call. = FALSE
      )
    }

    data_pp[, deltatime := as.numeric(tij)]
    data_pp <- data_pp[deltatime > 0]

    if ("node" %in% names(data_pp)) {
      data.table::setorderv(data_pp, c(".psl_pred_id", "node"))
    } else {
      data.table::setorderv(data_pp, ".psl_pred_id")
    }

    n_expanded <- nrow(data_pp)

    dt_pred <- vector("list", n_crisks)

    for (k in seq_len(n_crisks)) {

      pred_k <- object$model$private_predictor(
        model = object$learner_fit[[k]],
        newdata = endpoint_data,
        grid_nodes = grid_nodes
      )

      if (length(pred_k) != n_expanded) {
        stop(
          sprintf(
            "Prediction length mismatch for cause %s: got %s, expected %s.",
            k,
            length(pred_k),
            n_expanded
          ),
          call. = FALSE
        )
      }

      dt_pred[[k]] <- as.numeric(pred_k)
    }

    data_pp[, (pwch_cols) := dt_pred]

    haz <- as.matrix(data_pp[, .SD, .SDcols = pwch_cols])

    data_pp[
      ,
      survival_function := pch_survival(
        id = .psl_pred_id,
        dt = deltatime,
        haz = haz
      )
    ]

    data_pp[
      ,
      absolute_risk := pch_absolute_risk(
        .psl_pred_id,
        deltatime,
        haz,
        cause_idx = cause
      )
    ]

    last_pp <- data_pp[
      ,
      .SD[.N],
      by = .psl_pred_id
    ]

    pred_cols <- c(
      ".psl_pred_id",
      pwch_cols,
      "survival_function",
      "absolute_risk"
    )

    valid_out <- merge(
      valid_combo,
      last_pp[, ..pred_cols],
      by = ".psl_pred_id",
      all.x = TRUE,
      sort = FALSE
    )
  }


  ## ------------------------------------------------------------
  ## Assemble output
  ## ------------------------------------------------------------

  out <- data.table::rbindlist(
    list(zero_out, valid_out, future_out),
    use.names = TRUE,
    fill = TRUE
  )

  data.table::setorderv(
    out,
    c(".psl_internal_ix", ".psl_time_ix")
  )

  data.table::setnames(
    out,
    old = ".psl_internal_ix",
    new = id_col
  )

  private_cols <- c(".psl_pred_id", ".psl_time_ix")

  out[, (intersect(private_cols, names(out))) := NULL]

  return_cols <- unique(
    c(
      id_col,
      setdiff(
        names(newdata_dt),
        c(id_col, status_col, event_time_col, ".psl_internal_ix")
      ),
      event_time_col,
      pwch_cols,
      "survival_function",
      "absolute_risk"
    )
  )

  return_cols <- intersect(return_cols, names(out))

  out <- out[, ..return_cols]

  out[]
}
