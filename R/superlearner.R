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

  grid_nodes <- grid_nodes[grid_nodes <= max(data[[event_time]])]


  ## ------------------------------------------------------------
  ## Check learner library structure
  ## ------------------------------------------------------------

  is_learner <- function(x) {
    any(grepl("Learner_", class(x), fixed = TRUE))
  }

  is_learner_library <- function(x) {
    is.list(x) &&
      length(x) > 0L &&
      all(vapply(x, is_learner, logical(1L)))
  }

  if (!is.list(learners)) {
    stop(
      "'learners' must be either a single learner library or a list of learner libraries.",
      call. = FALSE
    )
  }

  if (is_learner_library(learners)) {

    if (n_crisks > 1L) {
      message(
        "A single learner library was supplied; the same learner library will be used for all competing causes."
      )
    }

    library_per_risk <- replicate(n_crisks, learners, simplify = FALSE)

  } else if (
    length(learners) == n_crisks &&
    all(vapply(learners, is_learner_library, logical(1L)))
  ) {

    library_per_risk <- learners

  } else {

    stop(
      paste0(
        "'learners' must be either:\n",
        "  (i) a single list of learners, which will be used for all competing causes, or\n",
        "  (ii) a list of ", n_crisks, " learner libraries, one for each competing cause."
      ),
      call. = FALSE
    )
  }


  ## ------------------------------------------------------------
  ## Name causes and learners
  ## ------------------------------------------------------------

  library_per_risk <- fill_missing_names(library_per_risk, prefix = "cause_")

  library_per_risk <- lapply(
    library_per_risk,
    fill_missing_names,
    prefix = "learner_"
  )


  grid_dt <- data.table(node = grid_nodes)
  grid_dt[, prev_node_psl := shift(node)]
  grid_dt[, (event_time) := node]

  dt <- grid_dt[data, on = event_time, roll = Inf]

  dt[, node_start := fifelse(get(event_time) == node, prev_node_psl, node)]


  dt <- dt[!is.na(node_start), c("node", "tij") := list(node_start, get(..event_time) - node_start)]


  setnames(dt, old = status, new = "deltaij")   ## deltaij is now different. It copies status and at each competing risk we need to evaluate deltaij := as.numeric(deltaij==competing_risk_j)

  dt[, folder := sample.int(nfold, nrow(dt), replace = TRUE)]

  ## ----------------------------------------------------------------
  ## Apply variable transformations in case the user provides them
  ## ----------------------------------------------------------------

  if (!is.null(variable_transformation)) {
    apply_transformations(dt, variable_transformation)
  }


  ## ------------------------------------------------------------
  ## Step 1: fit learners once on the full data, one list per cause
  ## ------------------------------------------------------------

  full_train_list <- vector("list", n_crisks)
  names(full_train_list) <- names(library_per_risk)

  for (jj in seq_len(n_crisks)) {
    full_train_list[[jj]] <- lapply(
      library_per_risk[[jj]],
      function(learner) {
        learner$private_fit(
          dt,
          cause = jj,
          grid_nodes = grid_nodes
        )
      }
    )

    names(full_train_list[[jj]]) <- names(library_per_risk[[jj]])
  }

  ## ------------------------------------------------------------
  ## Remove learners that fail on the full data
  ## ------------------------------------------------------------

  cause_labels <- names(library_per_risk)

  if (is.null(cause_labels) || anyNA(cause_labels) || any(!nzchar(cause_labels))) {
    cause_labels <- paste0("cause_", seq_along(library_per_risk))
  }

  failed_reason_table <- data.table::rbindlist(
    lapply(seq_along(full_train_list), function(jj) {

      fits_j <- full_train_list[[jj]]

      failed_j <- vapply(
        fits_j,
        is_failed_fit,
        logical(1L)
      )

      if (!any(failed_j)) {
        return(NULL)
      }

      data.table::data.table(
        cause_index = jj,
        cause = cause_labels[jj],
        learner_index = which(failed_j),
        learner = names(fits_j)[failed_j],
        reason = vapply(
          fits_j[failed_j],
          function(fit) {
            reason <- fit$reason

            if (is.null(reason) || length(reason) == 0L) {
              return(NA_character_)
            }

            as.character(reason)[1L]
          },
          character(1L)
        )
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )

  keep_by_cause <- lapply(
    full_train_list,
    function(fits_j) {
      !vapply(fits_j, is_failed_fit, logical(1L))
    }
  )

  n_retained_by_cause <- vapply(
    keep_by_cause,
    sum,
    integer(1L)
  )

  if (any(n_retained_by_cause == 0L)) {

    failed_causes <- cause_labels[n_retained_by_cause == 0L]

    reason_txt <- NULL

    if (nrow(failed_reason_table) > 0L) {
      reason_txt <- paste(
        apply(failed_reason_table, 1L, function(x) {
          sprintf(
            "%s, learner '%s': %s",
            x[["cause"]],
            x[["learner"]],
            x[["reason"]]
          )
        }),
        collapse = "\n"
      )
    }

    stop(
      paste(
        paste0(
          "All learners failed on the full data for at least one cause: ",
          paste(failed_causes, collapse = ", "),
          "."
        ),
        reason_txt,
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  if (nrow(failed_reason_table) > 0L) {

    warning_txt <- paste(
      apply(failed_reason_table, 1L, function(x) {
        sprintf(
          "%s, learner '%s' failed on the full data and was removed: %s",
          x[["cause"]],
          x[["learner"]],
          x[["reason"]]
        )
      }),
      collapse = "\n"
    )

    warning(warning_txt, call. = FALSE)

    for (jj in seq_along(library_per_risk)) {
      keep_j <- keep_by_cause[[jj]]

      library_per_risk[[jj]] <- library_per_risk[[jj]][keep_j]
      full_train_list[[jj]] <- full_train_list[[jj]][keep_j]
    }
  }

  learners_labels_per_risk <- lapply(library_per_risk, names)

  ## ------------------------------------------------------------
  ## If every cause has only one retained learner, return direct fit
  ## ------------------------------------------------------------

  n_learners_by_cause <- lengths(library_per_risk)

  if (all(n_learners_by_cause == 1L)) {

    learners_labels_per_risk <- lapply(library_per_risk, names)

    n_dropped_full <- if (
      exists("failed_reason_table") &&
      is.data.frame(failed_reason_table)
    ) {
      nrow(failed_reason_table)
    } else {
      0L
    }

    reason_txt <- NULL

    if (
      exists("failed_reason_table") &&
      is.data.frame(failed_reason_table) &&
      nrow(failed_reason_table) > 0L
    ) {
      reason_txt <- paste(
        sprintf(
          "%s, learner '%s': %s",
          failed_reason_table[["cause"]],
          failed_reason_table[["learner"]],
          failed_reason_table[["reason"]]
        ),
        collapse = "\n"
      )
    }

    msg <- paste0(
      "Only one usable base learner remains for every cause after full-data screening ",
      "(",
      n_dropped_full,
      " learner(s) dropped). Fitting each retained learner directly; ",
      "no cross-validation stacking is performed."
    )

    if (!is.null(reason_txt)) {
      msg <- paste(msg, reason_txt, sep = "\n")
    }

    message(msg)

    one_learner_out <- vector("list", n_crisks)
    names(one_learner_out) <- names(library_per_risk)

    for (cause_ix in seq_len(n_crisks)) {
      one_learner_out[[cause_ix]] <- list(
        learners_fit = full_train_list[[cause_ix]][[1L]],
        meta_learner_fit = NULL
      )
    }

    out <- list(
      learners = library_per_risk,
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
        learners_labels = learners_labels_per_risk,
        variable_transformation = variable_transformation
      )
    )

    class(out) <- "poisson_superlearner"

    return(out)
  }

  ## ------------------------------------------------------------
  ## Step 2: V-fold CV (only on retained learners)
  ## ------------------------------------------------------------

  # Create the level one data
  meta_learner_fits <- vector("list", n_crisks)
  cross_validation_deviance <- vector("list", n_crisks)


  ## ------------------------------------------------------------
  ## Precompute fold row indices once.
  ## ------------------------------------------------------------

  all_rows <- seq_len(nrow(dt))

  fold_rows <- vector("list", nfold)

  for (ix in seq_len(nfold)) {
    fold_rows[[ix]] <- which(dt[["folder"]] == ix)
  }


  for (jj in seq_len(n_crisks)) {

    risk_label_j <- names(library_per_risk)[jj]

    if (
      is.null(risk_label_j) ||
      is.na(risk_label_j) ||
      !nzchar(risk_label_j)
    ) {
      risk_label_j <- as.character(jj)
    }

    if (length(library_per_risk[[jj]]) == 1L) {
      meta_learner_fits[[jj]] <- NULL
      cross_validation_deviance[[jj]] <- NULL
      next
    }

    n_learners <- length(library_per_risk[[jj]])
    z_cols <- paste0("Z", seq_len(n_learners))

    fold_parts <- vector("list", nfold)

    ## Track learners that produce at least one non-NA prediction.
    ## This avoids scanning all Z columns after rbindlist().
    has_non_na <- logical(n_learners)


    for (ix in seq_len(nfold)) {

      valid_idx <- fold_rows[[ix]]

      valid_data <- dt[valid_idx]

      if (length(valid_idx)) {
        train_data <- dt[-valid_idx]
      } else {
        train_data <- dt[all_rows]
      }

      val_level <- make_validation_skeleton(
        valid_data = valid_data,
        grid_nodes = grid_nodes,
        cause = jj,
        node_col = "node",
        tij_col = "tij",
        event_col = "deltaij",

        ## If id/folder are needed later, change this back to:
        ## keep_cols = c(id, "folder")
        keep_cols = character(0L)
      )

      n_expanded <- nrow(val_level)

      train_list_jj <- lapply(
        library_per_risk[[jj]],
        function(ll) {
          ll$private_fit(
            train_data,
            cause = jj,
            grid_nodes = grid_nodes
          )
        }
      )

      preds <- vector("list", n_learners)

      for (mm in seq_len(n_learners)) {

        z <- library_per_risk[[jj]][[mm]]$private_predictor(
          model = train_list_jj[[mm]],
          newdata = valid_data,
          grid_nodes = grid_nodes
        )

        if (length(z) != n_expanded) {
          stop(
            sprintf(
              "Prediction length mismatch for cause %s, fold %s, learner %s: got %s, expected %s.",
              jj, ix, mm, length(z), n_expanded
            )
          )
        }

        z <- as.numeric(z)

        has_non_na[mm] <- has_non_na[mm] || any(!is.na(z))

        preds[[mm]] <- z
      }

      ## Assign all learner prediction columns in one data.table operation.
      val_level[, (z_cols) := preds]

      val_level[, fold := ix]

      fold_parts[[ix]] <- val_level

      rm(
        valid_idx,
        train_data,
        valid_data,
        val_level,
        train_list_jj,
        preds
      )
    }


    level_one_data <- data.table::rbindlist(
      fold_parts,
      use.names = TRUE,
      fill = FALSE
    )

    rm(fold_parts)


    ## ------------------------------------------------------------
    ## After CV, remove only learners that are ALL-NA for this cause.
    ## Partial fold failures are allowed.
    ## ------------------------------------------------------------

    learner_labels_j <- names(library_per_risk[[jj]])

    if (
      is.null(learner_labels_j) ||
      anyNA(learner_labels_j) ||
      any(!nzchar(learner_labels_j))
    ) {
      learner_labels_j <- paste0("learner_", seq_along(library_per_risk[[jj]]))
    }

    all_na_z <- !has_non_na

    if (all(all_na_z)) {

      failed_txt <- paste(
        sprintf(
          "cause %s, learner '%s': all cross-validated predictions are NA",
          risk_label_j,
          learner_labels_j
        ),
        collapse = "\n"
      )

      stop(
        paste(
          sprintf(
            "All learners failed during cross-validation for cause %s.",
            risk_label_j
          ),
          failed_txt,
          sep = "\n"
        ),
        call. = FALSE
      )
    }

    if (any(all_na_z)) {

      dropped_txt <- paste(
        sprintf(
          "cause %s, learner '%s': all cross-validated predictions are NA",
          risk_label_j,
          learner_labels_j[all_na_z]
        ),
        collapse = "\n"
      )

      warning(
        paste(
          "Removing learner(s) that produced only NA cross-validated predictions:",
          dropped_txt,
          sep = "\n"
        ),
        call. = FALSE
      )

      keep_ix <- which(!all_na_z)

      failed_z_cols <- z_cols[all_na_z]
      retained_old_z_cols <- z_cols[keep_ix]

      ## Drop failed prediction columns first, so renaming cannot collide.
      level_one_data[, (failed_z_cols) := NULL]

      ## Prune the cause-specific learner library and full-data fits.
      library_per_risk[[jj]] <- library_per_risk[[jj]][keep_ix]
      full_train_list[[jj]] <- full_train_list[[jj]][keep_ix]

      ## Re-index Z columns contiguously: Z1, Z2, ...
      new_z_cols <- paste0("Z", seq_along(retained_old_z_cols))

      data.table::setnames(
        level_one_data,
        old = retained_old_z_cols,
        new = new_z_cols
      )

      z_cols <- new_z_cols
      n_learners <- length(z_cols)

      learner_labels_j <- names(library_per_risk[[jj]])

      if (
        is.null(learner_labels_j) ||
        anyNA(learner_labels_j) ||
        any(!nzchar(learner_labels_j))
      ) {
        learner_labels_j <- paste0("learner_", seq_along(library_per_risk[[jj]]))
      }
    }


    ## ------------------------------------------------------------
    ## CV deviance on retained learners plus discrete superlearner
    ## ------------------------------------------------------------

    cv_dev_by_fold <- poisson_deviance_by_folder_hazard_cols(
      data = level_one_data,
      hazard_cols = z_cols,
      tij = as.numeric(level_one_data[["tij"]]),
      delta = as.integer(level_one_data[["deltaij"]]),
      fold = as.integer(level_one_data[["fold"]]),
      nfold = nfold
    )

    cv_dev <- colSums(cv_dev_by_fold)

    cv_dev[is.na(cv_dev) | !is.finite(cv_dev)] <- Inf
    names(cv_dev) <- names(library_per_risk[[jj]])

    cross_validation_deviance[[jj]] <- cv_dev

    ## -------------------------------------------------------------------
    ## Keep only relevant columns and transform to the logarithmic scale
    ## -------------------------------------------------------------------

    keep_cols <- c("tij", "deltaij", z_cols)
    drop_cols <- setdiff(names(level_one_data), keep_cols)

    if (length(drop_cols)) {
      level_one_data[, (drop_cols) := NULL]
    }

    level_one_data <- level_one_data[complete.cases(level_one_data)]

    eps <- 1e-15

    level_one_data[
      ,
      (z_cols) := lapply(
        .SD,
        function(x) log(pmax(as.numeric(x), eps))
      ),
      .SDcols = z_cols
    ]


    ## ------------------------------------------------------------
    ## If CV pruning leaves one learner for this cause, skip meta-learning
    ## ------------------------------------------------------------

    if (length(z_cols) == 1L) {
      meta_learner_fits[[jj]] <- NULL
      rm(level_one_data)
      next
    }


    ## ------------------------------------------------------------
    ## Meta-learning
    ## ------------------------------------------------------------

    x_meta <- as.matrix(level_one_data[, ..z_cols])
    storage.mode(x_meta) <- "double"

    meta_learner_fits[[jj]] <- glmnet::glmnet(
      x = x_meta,
      y = as.numeric(level_one_data[["deltaij"]]),
      offset = log(level_one_data[["tij"]]),
      family = "poisson",
      intercept = FALSE,
      lambda = 0
    )

    rm(level_one_data, x_meta)

    ## Optional, useful if this loop is memory-heavy.
    ## Avoid doing this inside each fold.
    gc()
  }

  #### OLD
  # for (jj in seq_len(n_crisks)) {
  #
  #
  #   if (length(library_per_risk[[jj]]) == 1L) {
  #     meta_learner_fits[[jj]] <- NULL
  #     cross_validation_deviance[[jj]] <- NULL
  #     next
  #   }
  #
  #   n_learners <- length(library_per_risk[[jj]])
  #   z_cols <- paste0("Z", seq_len(n_learners))
  #
  #   fold_parts <- vector("list", nfold)
  #
  #   for (ix in seq_len(nfold)) {
  #
  #     train_data <- dt[folder != ix]
  #     valid_data <- dt[folder == ix]
  #
  #     val_level <- make_validation_skeleton(
  #       valid_data = valid_data,
  #       grid_nodes = grid_nodes,
  #       cause = jj,
  #       node_col = "node",
  #       tij_col = "tij",
  #       event_col = "deltaij",
  #       keep_cols = c(id, "folder")
  #     )
  #
  #     n_expanded <- nrow(val_level)
  #
  #     train_list_jj <- lapply(
  #       library_per_risk[[jj]],
  #       function(ll) {
  #         ll$private_fit(
  #           train_data,
  #           cause = jj,
  #           grid_nodes = grid_nodes
  #         )
  #       }
  #     )
  #
  #     for (mm in seq_len(n_learners)) {
  #
  #       z <- library_per_risk[[jj]][[mm]]$private_predictor(
  #         model = train_list_jj[[mm]],
  #         newdata = valid_data,
  #         grid_nodes = grid_nodes
  #       )
  #
  #       if (length(z) != n_expanded) {
  #         stop(
  #           sprintf(
  #             "Prediction length mismatch for cause %s, fold %s, learner %s: got %s, expected %s.",
  #             jj, ix, mm, length(z), n_expanded
  #           )
  #         )
  #       }
  #
  #       val_level[, (z_cols[mm]) := as.numeric(z)]
  #     }
  #
  #
  #     val_level[, fold := ix]
  #
  #     fold_parts[[ix]] <- val_level
  #   }
  #
  #
  #   level_one_data <- data.table::rbindlist(
  #     fold_parts,
  #     use.names = TRUE,
  #     fill = FALSE
  #   )
  #
  #
  #
  #
  #   ## ------------------------------------------------------------
  #   ## After CV, remove only learners that are ALL-NA for this cause.
  #   ## Partial fold failures are allowed.
  #   ## ------------------------------------------------------------
  #
  #   learner_labels_j <- names(library_per_risk[[jj]])
  #
  #   if (
  #     is.null(learner_labels_j) ||
  #     anyNA(learner_labels_j) ||
  #     any(!nzchar(learner_labels_j))
  #   ) {
  #     learner_labels_j <- paste0("learner_", seq_along(library_per_risk[[jj]]))
  #   }
  #
  #   all_na_z <- vapply(
  #     z_cols,
  #     function(zc) all(is.na(level_one_data[[zc]])),
  #     logical(1L)
  #   )
  #
  #   if (all(all_na_z)) {
  #
  #     failed_txt <- paste(
  #       sprintf(
  #         "cause %s, learner '%s': all cross-validated predictions are NA",
  #         names(library_per_risk)[jj],
  #         learner_labels_j
  #       ),
  #       collapse = "\n"
  #     )
  #
  #     stop(
  #       paste(
  #         sprintf(
  #           "All learners failed during cross-validation for cause %s.",
  #           names(library_per_risk)[jj]
  #         ),
  #         failed_txt,
  #         sep = "\n"
  #       ),
  #       call. = FALSE
  #     )
  #   }
  #
  #   if (any(all_na_z)) {
  #
  #     dropped_txt <- paste(
  #       sprintf(
  #         "cause %s, learner '%s': all cross-validated predictions are NA",
  #         names(library_per_risk)[jj],
  #         learner_labels_j[all_na_z]
  #       ),
  #       collapse = "\n"
  #     )
  #
  #     warning(
  #       paste(
  #         "Removing learner(s) that produced only NA cross-validated predictions:",
  #         dropped_txt,
  #         sep = "\n"
  #       ),
  #       call. = FALSE
  #     )
  #
  #     keep_ix <- which(!all_na_z)
  #
  #     failed_z_cols <- z_cols[all_na_z]
  #     retained_old_z_cols <- z_cols[keep_ix]
  #
  #     ## Drop failed prediction columns first, so renaming cannot collide.
  #     level_one_data[, (failed_z_cols) := NULL]
  #
  #     ## Prune the cause-specific learner library and full-data fits.
  #     library_per_risk[[jj]] <- library_per_risk[[jj]][keep_ix]
  #     full_train_list[[jj]] <- full_train_list[[jj]][keep_ix]
  #
  #     ## Re-index Z columns contiguously: Z1, Z2, ...
  #     ## This avoids later prediction/meta-learner mismatches.
  #     new_z_cols <- paste0("Z", seq_along(retained_old_z_cols))
  #
  #     data.table::setnames(
  #       level_one_data,
  #       old = retained_old_z_cols,
  #       new = new_z_cols
  #     )
  #
  #     z_cols <- new_z_cols
  #     n_learners <- length(z_cols)
  #     learner_labels_j <- names(library_per_risk[[jj]])
  #   }
  #
  #   ## ------------------------------------------------------------
  #   ## CV deviance on retained learners plus discrete superlearner
  #   ## ------------------------------------------------------------
  #
  #   cv_dev_by_fold <- poisson_deviance_by_folder_hazard_cols(
  #     data = level_one_data,
  #     hazard_cols = z_cols,
  #     tij = as.numeric(level_one_data[["tij"]]),
  #     delta = as.integer(level_one_data[["deltaij"]]),
  #     fold = as.integer(level_one_data[["fold"]]),
  #     nfold = nfold
  #   )
  #
  #   cv_dev <- colSums(cv_dev_by_fold)
  #
  #   cv_dev[is.na(cv_dev) | !is.finite(cv_dev)] <- Inf
  #   names(cv_dev) <- names(library_per_risk[[jj]])
  #
  #   cross_validation_deviance[[jj]] <- cv_dev
  #
  #   discrete_ix <- which.min(cv_dev)
  #   discrete_col <- z_cols[discrete_ix]
  #
  #   ## -------------------------------------------------------------------
  #   ## Keep only relevant columns and transform to the logarithmic scale
  #   ## -------------------------------------------------------------------
  #
  #   level_one_data <- level_one_data[
  #     ,
  #     .SD,
  #     .SDcols = c("tij", "deltaij", z_cols)
  #   ]
  #
  #   level_one_data <- level_one_data[complete.cases(level_one_data), ]
  #
  #   eps <- 1e-15
  #
  #   level_one_data[
  #     ,
  #     (z_cols) := lapply(.SD, function(x) log(pmax(as.numeric(x), eps))),
  #     .SDcols = z_cols
  #   ]
  #
  #
  #   ## ------------------------------------------------------------
  #   ## If CV pruning leaves one learner for this cause, skip meta-learning
  #   ## ------------------------------------------------------------
  #
  #   if (length(z_cols) == 1L) {
  #     meta_learner_fits[[jj]] <- NULL
  #     next
  #   }
  #
  #
  #   ## ------------------------------------------------------------
  #   ## Meta-learning
  #   ## ------------------------------------------------------------
  #
  #   meta_learner_f <- create_formula_glmnet(
  #     covariates = z_cols,
  #     add_nodes = FALSE
  #   )
  #
  #   meta_learner_fits[[jj]] <- glmnet::glmnet(
  #     x = Matrix::sparse.model.matrix(
  #       formula(meta_learner_f),
  #       level_one_data,
  #       contrasts.arg = NULL
  #     )[, -1, drop = FALSE],
  #     y = as.numeric(level_one_data[["deltaij"]]),
  #     offset = log(level_one_data[["tij"]]),
  #     family = "poisson",
  #     intercept = FALSE,
  #     lambda = 0
  #   )
  #
  #
  # }


  ## ------------------------------------------------------------
  ## Step 3: build final output
  ## ------------------------------------------------------------

  names(meta_learner_fits) <- names(library_per_risk)
  names(cross_validation_deviance) <- names(library_per_risk)
  names(full_train_list) <- names(library_per_risk)

  learners_labels_per_risk <- lapply(library_per_risk, names)

  cross_validation_deviance_dt <- data.table::rbindlist(
    lapply(seq_along(cross_validation_deviance), function(jj) {

      cv_dev_j <- cross_validation_deviance[[jj]]

      if (is.null(cv_dev_j)) {
        return(NULL)
      }

      learner_labels_j <- names(cv_dev_j)

      if (
        is.null(learner_labels_j) ||
        anyNA(learner_labels_j) ||
        any(!nzchar(learner_labels_j))
      ) {
        learner_labels_j <- names(library_per_risk[[jj]])
      }

      data.table::data.table(
        cause_index = jj,
        cause = names(library_per_risk)[jj],
        learner_index = seq_along(cv_dev_j),
        learner = learner_labels_j,
        deviance = as.numeric(cv_dev_j)
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )

  superlearner_out <- lapply(seq_len(n_crisks), function(jj) {

    learners_fit_j <- full_train_list[[jj]]

    if (length(learners_fit_j) == 1L) {
      learners_fit_j <- learners_fit_j[[1L]]
    }

    list(
      learners_fit = learners_fit_j,
      meta_learner_fit = meta_learner_fits[[jj]]
    )
  })

  names(superlearner_out) <- names(library_per_risk)

  any_meta_learning <- any(
    !vapply(meta_learner_fits, is.null, logical(1L))
  )

  metalearner <- if (any_meta_learning) {
    list(
      engine = "glmnet::glmnet",
      family = "poisson",
      intercept = FALSE,
      lambda = 0,
      add_nodes = FALSE,
      scale = "log_hazard"
    )
  } else {
    NULL
  }

  if (!any_meta_learning) {
    message(
      "Only one usable base learner remains for every cause after cross-validation screening. No meta-learner was fitted."
    )
  }

  out <- list(
    learners = library_per_risk,
    metalearner = metalearner,
    superlearner = superlearner_out,
    cross_validation_deviance = cross_validation_deviance_dt,
    data_info = list(
      id = id,
      status = status,
      event_time = event_time,
      nodes = sort(unique(grid_nodes)),
      nfold = nfold,
      maximum_followup = maximum_followup,
      n_crisks = n_crisks,
      learners_labels = learners_labels_per_risk,
      variable_transformation = variable_transformation
    )
  )

  class(out) <- "poisson_superlearner"

  out

  }
