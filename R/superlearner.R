#' Fit a Poisson Super Learner ensemble
#'
#' Fits an ensemble of cause-specific piecewise-constant hazard models using a
#' long-format Poisson representation and combines them through a meta-learner
#' (stacking).
#'
#' Internally, the function:
#' \enumerate{
#'   \item builds a time grid (`nodes`) and converts the subject-level data to a
#'   long Poisson format;
#'   \item fits each base learner once on the full long data for each cause;
#'   \item removes learners that already fail on the full data;
#'   \item uses `nfold` cross-validation to obtain out-of-sample base-learner
#'   predictions (`Z1`, `Z2`, ...) for stacking;
#'   \item removes learners whose cross-validated prediction column is entirely
#'   missing for at least one cause;
#'   \item fits a cause-specific meta-learner on the retained stacked predictions.
#' }
#'
#' @param data `data.frame`. Subject-level input data, one row per subject.
#' @param id `character(1)`. Name of the subject identifier column. If missing,
#'   an `id` column is created automatically.
#' @param status `character(1)`. Name of the event-status column. It must be coded
#'   with `0` for censoring and `1, 2, ..., K` for event types. If there is no
#'   `0` in `status`, the data are treated as uncensored.
#' @param event_time `character(1)`. Name of the event or censoring time column.
#' @param learners `list`. List of initialized learner reference-class objects,
#'   for example [Learner_glmnet()], [Learner_hal()], or [Learner_gam()]. If
#'   unnamed, learners are named `"learner_1"`, `"learner_2"`, and so on. Each
#'   learner must implement `$private_fit(dt_long)` and
#'   `$private_predictor(model, newdata)`.
#' @param number_of_nodes `numeric(1)` or `NULL`. If not `NULL`, constructs a
#'   quantile-based node grid with `number_of_nodes + 1` cut points. Ignored when
#'   `nodes` is supplied.
#' @param nodes `numeric` or `NULL`. Explicit time-node grid. If supplied,
#'   `number_of_nodes` is ignored. `0` is added if missing, and nodes larger than
#'   `max(event_time)` are dropped.
#' @param variable_transformation Optional transformation specification passed to
#'   `apply_transformations()` on the internally created long-format data.
#' @param nfold `numeric(1)`. Number of folds for cross-validation stacking.
#' @param ... Additional arguments currently ignored.
#'
#' @return An object of class `poisson_superlearner`, stored as a named `list`
#'   with the following components:
#'
#'   `learners`:
#'   the retained base learner objects.
#'
#'   `metalearner`:
#'   the meta-learner object used for stacking. If no stacking is performed because
#'   only one learner remains, `metalearner` is `NULL`.
#'
#'   `superlearner`:
#'   a `list` of length `data_info$n_crisks`, one entry per cause. For cause `k`,
#'   `superlearner[[k]]` is a `list` with two elements:
#'   \itemize{
#'     \item `learners_fit`: the fitted base learner object or objects for cause `k`.
#'     If more than one learner is retained, this is a `list` with one fitted
#'     object per retained learner. If only one learner remains, this is the
#'     single fitted learner object itself.
#'     \item `meta_learner_fit`: the fitted cause-specific meta-learner for cause
#'     `k`. If no stacking is performed, this is `NULL`.
#'   }
#'
#'   `cross_validation_deviance`:
#'   a `data.table` with columns `learner` and `deviance`, giving the mean
#'   cross-validated Poisson deviance for each retained base learner. This
#'   component is present when cross-validated model comparison is available.
#'
#'   `data_info`:
#'   a `list` of bookkeeping information used for prediction and interpretation,
#'   containing:
#'   \itemize{
#'     \item `id`: identifier column name used.
#'     \item `status`: status column name used.
#'     \item `event_time`: event-time column name used.
#'     \item `nodes`: numeric vector of node cut points used for the piecewise grid.
#'     \item `nfold`: number of folds used for stacking.
#'     \item `maximum_followup`: maximum observed follow-up time.
#'     \item `n_crisks`: number of event types detected.
#'     \item `learners_labels`: character vector of retained learner labels.
#'     \item `variable_transformation`: the transformation specification passed in
#'     `variable_transformation`, or `NULL`.
#'   }
#'
#' @details
#' If all learners fail on the full data, the function stops with an error.
#' If only one learner remains after the full-data screening step or after the
#' cross-validation screening step, no meta-learner is fit. In that case,
#' `metalearner` is `NULL`, each `superlearner[[k]]$meta_learner_fit` is `NULL`,
#' and prediction is based directly on the stored fitted base learner.Numeric
#' learner positions always refer to the learners actually retained in the
#' fitted object.
#'
#' @examples
#' data <- simulateStenoT1(50, competing_risks = TRUE)
#'
#' learners <- list(
#'   glm = Learner_glmnet(
#'     covariates = c("sex", "value_LDL"),
#'     lambda = 0,
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
#'   data = data,
#'   id = "id",
#'   status = "status_cvd",
#'   event_time = "time_cvd",
#'   learners = learners,
#'   number_of_nodes = 3,
#'   nfold = 2
#' )
#'
#' @export
Superlearner <- function(data,
                         id = "id",
                         status = "status",
                         event_time = NULL,
                         learners,
                         number_of_nodes = NULL,
                         nodes = NULL,
                         variable_transformation = NULL,
                         nfold = 3,
                         ...) {

  if (!(id %in% names(data))) {
    data[["id"]] <- 1:NROW(data)
    id <- "id"
  }

  # if (is.null(names(learners))) {
  #   names(learners) <- paste0("learner_", seq_along(learners))
  # }


  # learners_labels <- names(learners)

  maximum_followup <- max(data[[event_time]])
  n <- length(unique(data[[id]]))

  if (!(0 %in% data[[status]])) {
    warning(
      paste0(
        "There is no value of ",
        status,
        " equal to zero. We will consider the data uncensored."
      )
    )
    n_crisks <- length(unique(data[[status]]))
    uncensored_01 <- TRUE
  } else {
    n_crisks <- length(unique(data[[status]])) - 1L
    uncensored_01 <- FALSE
  }

  if (!is.null(number_of_nodes)) {
    grid_nodes <- quantile(
      data[[event_time]],
      probs = seq(0, 1, length.out = as.integer(number_of_nodes) + 1L),
      type = 1,
      names = FALSE
    )
  } else if (is.null(nodes)) {
    grid_nodes <- sort(unique(data[[event_time]]))
  } else {
    grid_nodes <- sort(nodes)
  }

  if (!(0 %in% grid_nodes)) {
    grid_nodes <- c(0, grid_nodes)
  }

  grid_nodes <- grid_nodes[grid_nodes <= maximum_followup]

################################
  is_psl_learner <- function(x) {
    is.environment(x) &&
      is.function(x$private_fit) &&
      is.function(x$private_predictor)
  }

  is_flat_library <- is.list(learners) &&
    length(learners) > 0L &&
    all(vapply(learners, is_psl_learner, logical(1)))

  is_cause_specific_library <- is.list(learners) &&
    length(learners) == n_crisks &&
    all(vapply(learners, function(x) {
      is.list(x) &&
        length(x) > 0L &&
        all(vapply(x, is_psl_learner, logical(1)))
    }, logical(1)))

  if (is_flat_library) {
    if (is.null(names(learners))) {
      names(learners) <- paste0("learner_", seq_along(learners))
    }

    learners_by_cause <- replicate(n_crisks, learners, simplify = FALSE)

  } else if (is_cause_specific_library) {
    learners_by_cause <- learners

    for (k in seq_len(n_crisks)) {
      if (is.null(names(learners_by_cause[[k]]))) {
        names(learners_by_cause[[k]]) <- paste0("learner_", seq_along(learners_by_cause[[k]]))
      }
    }

  } else {
    stop(
      "`learners` must either be a list of learner objects, ",
      "or a list of length n_crisks where each element is a list of learner objects.",
      call. = FALSE
    )
  }

  learners_labels_by_cause <- lapply(learners_by_cause, names)
  z_covariates_by_cause <- lapply(
    learners_by_cause,
    function(x) paste0("Z", seq_along(x))
  )

  ################################

  dt <- data_pre_processing(
    data = data,
    id = id,
    status = status,
    nodes = grid_nodes,
    event_time = event_time,
    uncensored_01 = uncensored_01
  )

  if (!is.null(variable_transformation)) {
    apply_transformations(dt, variable_transformation)
  }

  id_fold <- sample(
    1:nfold,
    n,
    replace = TRUE,
    prob = rep(1 / nfold, nfold)
  )

  dt_id <- data.table(folder = id_fold, id = unique(data[[id]]))
  dt <- merge(dt, dt_id, by = "id", all.x = TRUE)

  dt_by_cause <- split(dt, by = "k")

  # z_covariates <- paste0("Z", seq_along(learners))

  ## ------------------------------------------------------------
  ## Step 1: fit learners once on the full data, one list per cause
  ## ------------------------------------------------------------
  # full_train_list <- lapply(dt_by_cause, function(dt_k) {
  #   lapply(learners, function(f) f$private_fit(dt_k))
  # })
  full_train_list <- mapply(
    function(dt_k, learners_k) {
      lapply(learners_k, function(f) f$private_fit(dt_k))
    },
    dt_by_cause,
    learners_by_cause,
    SIMPLIFY = FALSE
  )
  ## ------------------------------------------------------------
  ## Save failure reasons before pruning, cause-specific version
  ## ------------------------------------------------------------

  failed_reason_table <- data.table::rbindlist(
    lapply(seq_along(full_train_list), function(k) {
      fits_k <- full_train_list[[k]]

      data.table::data.table(
        cause = k,
        learner = learners_labels_by_cause[[k]],
        z = z_covariates_by_cause[[k]],
        reason = vapply(
          fits_k,
          function(fit) {
            if (is_failed_fit(fit)) fit$reason else NA_character_
          },
          character(1)
        )
      )
    }),
    fill = TRUE
  )

  failed_reason_table <- failed_reason_table[
    !is.na(reason)
  ]

  reason_txt <- NULL
  if (nrow(failed_reason_table) > 0L) {
    reason_txt <- paste(
      failed_reason_table[
        ,
        sprintf(
          "cause %s, learner '%s': %s",
          cause,
          learner,
          reason
        )
      ],
      collapse = "\n"
    )
  }


  ## ------------------------------------------------------------
  ## Remove full-data failures within each cause only
  ## ------------------------------------------------------------

  for (k in seq_len(n_crisks)) {
    failed_k <- vapply(
      full_train_list[[k]],
      is_failed_fit,
      logical(1)
    )

    keep_k <- !failed_k

    learners_by_cause[[k]] <- learners_by_cause[[k]][keep_k]
    learners_labels_by_cause[[k]] <- learners_labels_by_cause[[k]][keep_k]
    full_train_list[[k]] <- full_train_list[[k]][keep_k]

    ## Re-index Z columns within cause after pruning.
    ## This is simpler and safer than carrying gaps like Z1, Z3.
    z_covariates_by_cause[[k]] <- paste0(
      "Z",
      seq_along(learners_by_cause[[k]])
    )
  }


  ## ------------------------------------------------------------
  ## Stop only if a cause has no usable learner
  ## ------------------------------------------------------------

  n_learners_by_cause <- lengths(learners_by_cause)
  empty_causes <- which(n_learners_by_cause == 0L)

  if (length(empty_causes) > 0L) {
    stop(
      paste(
        paste0(
          "All learners failed on the full data before cross-validation for cause(s): ",
          paste(empty_causes, collapse = ", "),
          "."
        ),
        reason_txt,
        sep = "\n"
      ),
      call. = FALSE
    )
  }


  ## ------------------------------------------------------------
  ## Flags for the next steps
  ## ------------------------------------------------------------

  needs_meta_by_cause <- n_learners_by_cause > 1L
  all_direct_by_cause <- all(n_learners_by_cause == 1L)
  any_direct_by_cause <- any(n_learners_by_cause == 1L)

  ## ------------------------------------------------------------
  ## If every cause has only one learner, return direct fit
  ## ------------------------------------------------------------

  if (all_direct_by_cause) {

    msg <- paste0(
      "Only one usable base learner remains for every cause after full-data screening. ",
      "Fitting cause-specific learners directly; no ensemble constructed."
    )

    if (!is.null(reason_txt)) {
      msg <- paste(msg, reason_txt, sep = "\n")
    }

    message(msg)

    one_learner_out <- vector("list", n_crisks)

    for (cause_ix in seq_len(n_crisks)) {
      one_learner_out[[cause_ix]] <- list(
        learners_fit = full_train_list[[cause_ix]][[1L]],
        meta_learner_fit = NULL
      )
    }

    out <- list(
      ## New source of truth
      learners_by_cause = learners_by_cause,
      learners_labels_by_cause = learners_labels_by_cause,
      z_covariates_by_cause = z_covariates_by_cause,

      ## Backward-compatible alias only when the original input was one flat library
      learners = if (is_flat_library) learners_by_cause[[1L]] else NULL,

      metalearner = NULL,
      superlearner = one_learner_out,

      data_info = list(
        id = id,
        status = status,
        event_time = event_time,
        nodes = sort(unique(grid_nodes)),
        nfold = nfold,
        maximum_followup = maximum_followup,
        n_crisks = n_crisks,

        ## New bookkeeping
        learners_labels_by_cause = learners_labels_by_cause,
        n_learners_by_cause = n_learners_by_cause,

        ## Backward-compatible alias only when meaningful
        learners_labels = if (is_flat_library) learners_labels_by_cause[[1L]] else NULL,

        variable_transformation = variable_transformation
      )
    )

    class(out) <- "poisson_superlearner"
    return(out)
  }


  ## ------------------------------------------------------------
  ## Step 2: V-fold CV only on retained learners
  ## ------------------------------------------------------------

  dt_by_cause_cv <- lapply(dt_by_cause, function(dt_k) {
    out <- data.table::copy(dt_k)
    out[, row_ix := .I]
    out
  })

  ## One OOF buffer per cause, with cause-specific Z columns
  oof_buffers <- mapply(
    function(dt_k, z_k) {
      psl_init_oof_buffer(
        dt_k = dt_k,
        z_covariates = z_k
      )
    },
    dt_k = dt_by_cause_cv,
    z_k = z_covariates_by_cause,
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )

  for (ix in seq_len(nfold)) {

    pseudo_observations <- mapply(
      function(dt_k, competing_risk, learners_k, z_k) {
        create_pseudo_observations(
          training_data = dt_k[folder != ix, ],
          validation_data = dt_k[folder == ix, ],
          competing_risk = competing_risk,
          learners = learners_k,
          z_covariates = z_k,
          ix = ix
        )
      },
      dt_k = dt_by_cause_cv,
      competing_risk = as.list(seq_len(n_crisks)),
      learners_k = learners_by_cause,
      z_k = z_covariates_by_cause,
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )

    for (k in seq_len(n_crisks)) {
      oof_buffers[[k]] <- psl_write_oof_buffer(
        oof_buffer = oof_buffers[[k]],
        oof_chunk = pseudo_observations[[k]]
      )
    }
  }

  ## ------------------------------------------------------------
  ## Step 3: CV deviance on retained learners, cause-specific version
  ## ------------------------------------------------------------

  dt_learners_by_cause <- vector("list", n_crisks)

  for (k in seq_len(n_crisks)) {

    z_k <- z_covariates_by_cause[[k]]

    loghaz_cols <- psl_oof_loghaz_dt(
      oof_buffer = oof_buffers[[k]],
      z_covariates = z_k
    )

    dev_k <- poisson_deviance_by_folder_cols(
      log_hazard_cols = loghaz_cols,
      tij = as.numeric(dt_by_cause_cv[[k]][["tij"]]),
      delta = as.integer(dt_by_cause_cv[[k]][["deltaij"]]),
      folder = as.integer(dt_by_cause_cv[[k]][["folder"]]),
      nfold = nfold
    )

    dev_mean_k <- colMeans(2.0 * dev_k)

    dt_learners_by_cause[[k]] <- data.table::data.table(
      cause = k,
      learner = learners_labels_by_cause[[k]],
      z = z_k,
      deviance = dev_mean_k
    )
  }

  dt_learners <- data.table::rbindlist(
    dt_learners_by_cause,
    use.names = TRUE,
    fill = TRUE
  )

  ## Keep only the identifying columns needed to reattach deviance after pruning.
  ## Do NOT keep old z here, because z may be re-indexed after pruning.
  cross_validation_deviance_before_cv_pruning <- data.table::copy(
    dt_learners[, .(cause, learner, deviance)]
  )

  ## Safety: learner labels must be unique within each cause
  duplicated_labels_by_cause <- lapply(seq_len(n_crisks), function(k) {
    learners_labels_by_cause[[k]][duplicated(learners_labels_by_cause[[k]])]
  })

  bad_causes <- which(lengths(duplicated_labels_by_cause) > 0L)

  if (length(bad_causes) > 0L) {
    stop(
      paste0(
        "Learner labels must be unique within each cause. Duplicated labels found in cause(s): ",
        paste(bad_causes, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }


  ## ------------------------------------------------------------
  ## Step 4: after CV, remove only learners that are ALL-NA
  ## within their own cause. Partial fold failures are allowed.
  ## ------------------------------------------------------------

  failed_cv_by_cause <- vector("list", n_crisks)

  for (k in seq_len(n_crisks)) {

    z_k <- z_covariates_by_cause[[k]]
    buf_k <- oof_buffers[[k]]

    failed_k <- colSums(!is.na(buf_k$Z)) == 0L
    failed_cv_by_cause[[k]] <- z_k[failed_k]

    if (any(failed_k)) {

      keep_ix <- which(!failed_k)

      learners_by_cause[[k]] <- learners_by_cause[[k]][keep_ix]
      learners_labels_by_cause[[k]] <- learners_labels_by_cause[[k]][keep_ix]
      full_train_list[[k]] <- full_train_list[[k]][keep_ix]

      buf_k$Z <- buf_k$Z[, keep_ix, drop = FALSE]

      ## Re-index Z columns after pruning inside this cause.
      ## This avoids gaps such as Z1, Z3.
      z_covariates_by_cause[[k]] <- paste0(
        "Z",
        seq_along(learners_by_cause[[k]])
      )

      colnames(buf_k$Z) <- z_covariates_by_cause[[k]]
      oof_buffers[[k]] <- buf_k
    }
  }


  ## ------------------------------------------------------------
  ## Stop if any cause has no usable learner left
  ## ------------------------------------------------------------

  n_learners_by_cause <- lengths(learners_by_cause)
  empty_causes <- which(n_learners_by_cause == 0L)

  if (length(empty_causes) > 0L) {
    stop(
      paste0(
        "All learners failed after cross-validation screening for cause(s): ",
        paste(empty_causes, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }


  ## ------------------------------------------------------------
  ## Rebuild the CV deviance table after CV pruning
  ## ------------------------------------------------------------

  dt_learners_by_cause <- lapply(seq_len(n_crisks), function(k) {
    data.table::data.table(
      cause = k,
      learner = learners_labels_by_cause[[k]],
      z = z_covariates_by_cause[[k]]
    )
  })

  dt_learners <- data.table::rbindlist(
    dt_learners_by_cause,
    use.names = TRUE,
    fill = TRUE
  )

  dt_learners <- merge(
    dt_learners,
    cross_validation_deviance_before_cv_pruning,
    by = c("cause", "learner"),
    all.x = TRUE,
    sort = FALSE
  )

  data.table::setorder(dt_learners, cause)


  ## ------------------------------------------------------------
  ## Step 5: meta-learning, cause-specific version
  ## ------------------------------------------------------------

  n_learners_by_cause <- lengths(learners_by_cause)
  needs_meta_by_cause <- n_learners_by_cause > 1L
  all_direct_by_cause <- all(n_learners_by_cause == 1L)

  ## Defensive: if this was not already handled after Step 4
  if (all_direct_by_cause) {

    message(
      "Only one usable base learner remains for every cause after cross-validation screening. ",
      "No ensemble constructed."
    )

    superlearner_out <- vector("list", n_crisks)

    for (cause_ix in seq_len(n_crisks)) {
      superlearner_out[[cause_ix]] <- list(
        learners_fit = full_train_list[[cause_ix]][[1L]],
        meta_learner_fit = NULL
      )
    }

    same_retained_library <- is_flat_library &&
      length(unique(vapply(
        learners_labels_by_cause,
        paste,
        character(1),
        collapse = "\r"
      ))) == 1L

    learners_alias <- if (same_retained_library) learners_by_cause[[1L]] else NULL
    learners_labels_alias <- if (same_retained_library) learners_labels_by_cause[[1L]] else NULL

    out <- list(
      learners_by_cause = learners_by_cause,
      learners_labels_by_cause = learners_labels_by_cause,
      z_covariates_by_cause = z_covariates_by_cause,

      ## Backward-compatible aliases only when meaningful
      learners = learners_alias,

      metalearner_by_cause = vector("list", n_crisks),
      metalearner = NULL,

      superlearner = superlearner_out,
      cross_validation_deviance = dt_learners,

      data_info = list(
        id = id,
        status = status,
        event_time = event_time,
        nodes = sort(unique(grid_nodes)),
        nfold = nfold,
        maximum_followup = maximum_followup,
        n_crisks = n_crisks,

        learners_labels_by_cause = learners_labels_by_cause,
        z_covariates_by_cause = z_covariates_by_cause,
        n_learners_by_cause = n_learners_by_cause,

        learners_labels = learners_labels_alias,

        variable_transformation = variable_transformation
      )
    )

    class(out) <- "poisson_superlearner"
    return(out)
  }


  ## ------------------------------------------------------------
  ## Fit one meta-learner per cause that still has >1 learner
  ## ------------------------------------------------------------

  meta_learner_by_cause <- vector("list", n_crisks)
  meta_learner_fits <- vector("list", n_crisks)

  for (k in seq_len(n_crisks)) {

    if (!needs_meta_by_cause[k]) {
      meta_learner_by_cause[[k]] <- NULL
      meta_learner_fits[[k]] <- NULL
      next
    }

    z_k <- z_covariates_by_cause[[k]]

    meta_learner_by_cause[[k]] <- Learner_glmnet(
      covariates = z_k,
      cross_validation = FALSE,
      intercept = FALSE,
      add_nodes = FALSE,
      penalise_nodes = TRUE,
      lambda = 0
    )

    meta_learner_fits[[k]] <- fit_meta_learner(
      dt = dt_by_cause[[k]],
      dt_z = oof_buffers[[k]],
      meta_learner = meta_learner_by_cause[[k]],
      z_covariates = z_k
    )
  }


  ## ------------------------------------------------------------
  ## Store cause-specific learner fits and meta-learner fits
  ## ------------------------------------------------------------

  superlearner_out <- vector("list", n_crisks)

  for (k in seq_len(n_crisks)) {

    if (n_learners_by_cause[k] == 1L) {
      superlearner_out[[k]] <- list(
        learners_fit = full_train_list[[k]][[1L]],
        meta_learner_fit = NULL
      )
    } else {
      superlearner_out[[k]] <- list(
        learners_fit = full_train_list[[k]],
        meta_learner_fit = meta_learner_fits[[k]]
      )
    }
  }


  ## ------------------------------------------------------------
  ## Backward-compatible aliases
  ## ------------------------------------------------------------

  same_retained_library <- is_flat_library &&
    length(unique(vapply(
      learners_labels_by_cause,
      paste,
      character(1),
      collapse = "\r"
    ))) == 1L

  learners_alias <- if (same_retained_library) learners_by_cause[[1L]] else NULL
  learners_labels_alias <- if (same_retained_library) learners_labels_by_cause[[1L]] else NULL

  ## Only meaningful if all causes that need meta-learning use the same Z layout.
  same_meta_layout <- any(needs_meta_by_cause) &&
    length(unique(vapply(
      z_covariates_by_cause[needs_meta_by_cause],
      paste,
      character(1),
      collapse = "\r"
    ))) == 1L

  metalearner_alias <- if (same_meta_layout) {
    meta_learner_by_cause[[which(needs_meta_by_cause)[1L]]]
  } else {
    NULL
  }


  ## ------------------------------------------------------------
  ## Final object
  ## ------------------------------------------------------------

  out <- list(
    ## New source of truth
    learners_by_cause = learners_by_cause,
    learners_labels_by_cause = learners_labels_by_cause,
    z_covariates_by_cause = z_covariates_by_cause,

    metalearner_by_cause = meta_learner_by_cause,

    ## Backward-compatible aliases only when meaningful
    learners = learners_alias,
    metalearner = metalearner_alias,

    superlearner = superlearner_out,
    cross_validation_deviance = dt_learners,

    data_info = list(
      id = id,
      status = status,
      event_time = event_time,
      nodes = sort(unique(grid_nodes)),
      nfold = nfold,
      maximum_followup = maximum_followup,
      n_crisks = n_crisks,

      learners_labels_by_cause = learners_labels_by_cause,
      z_covariates_by_cause = z_covariates_by_cause,
      n_learners_by_cause = n_learners_by_cause,
      needs_meta_by_cause = needs_meta_by_cause,

      ## Old-style label vector, only if still unambiguous
      learners_labels = learners_labels_alias,

      variable_transformation = variable_transformation
    )
  )

  class(out) <- "poisson_superlearner"
  out
}
