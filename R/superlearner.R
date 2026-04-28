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

  for (jj in seq_len(n_crisks)) {
    full_train_list[[jj]] <- lapply(
      learners[[jj]],
      function(learner) {
        learner$private_fit(
          dt,
          cause = jj,
          grid_nodes = grid_nodes
        )
      }
    )
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

  ##########################################################

  ## If pruning leaves one learner, return direct fit
  # if (length(learners) == 1L) {
  #
  #   msg <- paste0(
  #     "Only one usable base learner remains after full-data screening (",
  #     nrow(failed_reason_table),
  #     " learner(s) dropped). Fitting the learner directly; no ensemble constructed."
  #   )
  #
  #   if (!is.null(reason_txt)) {
  #     msg <- paste(msg, reason_txt, sep = "\n")
  #   }
  #
  #   message(msg)
  #
  #   one_learner_out <- vector("list", n_crisks)
  #
  #   for (cause_ix in seq_len(n_crisks)) {
  #     one_learner_out[[cause_ix]] <- list(
  #       learners_fit = full_train_list[[cause_ix]][[1L]],
  #       meta_learner_fit = NULL
  #     )
  #   }
  #
  #   out <- list(
  #     learners = learners,
  #     metalearner = NULL,
  #     superlearner = one_learner_out,
  #     data_info = list(
  #       id = id,
  #       status = status,
  #       event_time = event_time,
  #       # nodes = sort(unique(as.numeric(levels(dt$node)))),
  #       nodes = sort(unique(grid_nodes)),
  #       nfold = nfold,
  #       maximum_followup = maximum_followup,
  #       n_crisks = n_crisks,
  #       learners_labels = learners_labels,
  #       variable_transformation = variable_transformation
  #     )
  #   )
  #
  #   class(out) <- "poisson_superlearner"
  #   return(out)
  # }

  ## ------------------------------------------------------------
  ## Step 2: V-fold CV only on retained learners
  ## ------------------------------------------------------------

  # Create the level one data
  meta_learner_fits <- vector("list", n_crisks)
  cross_validation_deviance <- vector("list", n_crisks)

  for (jj in seq_len(n_crisks)) {

    n_learners <- length(library_per_risk[[jj]])
    z_cols <- paste0("Z", seq_len(n_learners))

    fold_parts <- vector("list", nfold)

    for (ix in seq_len(nfold)) {

      train_data <- dt[folder != ix]
      valid_data <- dt[folder == ix]

      val_level <- make_validation_skeleton(
        valid_data = valid_data,
        grid_nodes = grid_nodes,
        cause = jj,
        node_col = "node",
        tij_col = "tij",
        event_col = "deltaij",
        keep_cols = c(id, "folder")
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

        val_level[, (z_cols[mm]) := as.numeric(z)]
      }


      val_level[, fold := ix]

      fold_parts[[ix]] <- val_level
    }


    level_one_data <- data.table::rbindlist(
      fold_parts,
      use.names = TRUE,
      fill = FALSE
    )

    ## ------------------------------------------------------------
    ## CV deviance on retained learners plus discrete superlearner
    ## ------------------------------------------------------------

    all_na_z <- vapply(
      z_cols,
      function(zc) all(is.na(level_one_data[[zc]])),
      logical(1)
    )

    retained_z_cols <- z_cols[!all_na_z]

    if (!length(retained_z_cols)) {
      stop(sprintf("All learners failed for cause %s.", jj))
    }

    cv_dev_by_fold <- poisson_deviance_by_folder_hazard_cols(
      data = level_one_data,
      hazard_cols = retained_z_cols,
      tij = as.numeric(level_one_data[["tij"]]),
      delta = as.integer(level_one_data[["deltaij"]]),
      fold = as.integer(level_one_data[["fold"]]),
      nfold = nfold
    )

    cv_dev <- colSums(cv_dev_by_fold)

    cv_dev[is.na(cv_dev) | !is.finite(cv_dev)] <- Inf

    cross_validation_deviance[[jj]] <- cv_dev

    discrete_ix <- which.min(cv_dev)
    discrete_col <- retained_z_cols[discrete_ix]

    ## -------------------------------------------------------------------
    ## Keep only relevant columns and transform to the logarithmic scale
    ## -------------------------------------------------------------------


    level_one_data<-level_one_data[,.SD,.SDcols=c("tij","deltaij",z_cols)]

    eps <- 1e-15

    level_one_data[
      ,
      (z_cols) := lapply(.SD, function(x) log(pmax(as.numeric(x), eps))),
      .SDcols = z_cols
    ]

    ## ------------------------------------------------------------
    ## After CV, remove only learners that are ALL-NA
    ## within at least one cause. Partial fold failures are allowed.
    ## ------------------------------------------------------------


    ## ------------------------------------------------------------
    ## Change of logic: Meta-learning, without duplicating the meta-learner.
    ## ------------------------------------------------------------


    meta_learner_f <- create_formula_glmnet(covariates = z_cols, add_nodes =
                                              FALSE)


    meta_learner_fits[[jj]] <- glmnet(
      x = Matrix::sparse.model.matrix(formula(meta_learner_f), level_one_data, contrasts.arg = NULL)[, -1, drop = FALSE],
      y = as.numeric(level_one_data[["deltaij"]]),
      offset = log(level_one_data[["tij"]]),
      family = "poisson",
      intercept = FALSE,
      lambda = 0
    )


  }


browser()




    # for(jj in seq_along(n_crisks)){
    #
    #   ## create a structure with nrows == n expandend rows and
    #   ## ncols == length(library_per_risk[[jj]]) + 2 (deltaij expanded and tij expanded)
    #
    #   for (ix in seq_len(nfold)) {
    #   browser()
    #
    #   train_list_jj <- lapply(library_per_risk[[jj]], function(ll)
    #     ll$private_fit(dt[folder != ix, ], cause = jj, grid_nodes = grid_nodes))
    #
    #   val_list <- mapply(
    #     function(ll, model, newdata, grid_nodes) {
    #       ll$private_predictor(model = model,
    #                           newdata = newdata,
    #                           grid_nodes=grid_nodes)
    #     },
    #     library_per_risk[[jj]],
    #     train_list_jj,
    #     MoreArgs = list(newdata = dt[folder == ix, ],
    #                     grid_nodes=grid_nodes),
    #     SIMPLIFY = FALSE
    #   )
    #
    #   ### save into a structure with expanded dimensions
    #
    #   }
    #
    #   #### reduce the dimensions to nnodes x nlearners
    #
    # }






    # pseudo_observations <- mapply(
    #   function(dt_k, competing_risk) {
    #     create_pseudo_observations(
    #       training_data = dt_k[folder != ix, ],
    #       validation_data = dt_k[folder == ix, ],
    #       competing_risk = competing_risk,
    #       learners = learners,
    #       z_covariates = z_covariates,
    #       ix = ix
    #     )
    #   },
    #   dt_k = dt_by_cause_cv,
    #   competing_risk = as.list(seq_len(n_crisks)),
    #   SIMPLIFY = FALSE,
    #   USE.NAMES = FALSE
    # )
    #
    # for (k in seq_len(n_crisks)) {
    #   oof_buffers[[k]] <- psl_write_oof_buffer(
    #     oof_buffer = oof_buffers[[k]],
    #     oof_chunk = pseudo_observations[[k]]
    #   )
    # }
  # }

## ------------------------------------------------------------
## Step 3: CV deviance on retained learners
## ------------------------------------------------------------




# L <- length(z_covariates)
# dev_sum <- matrix(0.0, nrow = nfold, ncol = L)
#
# for (k in seq_len(n_crisks)) {
#   loghaz_cols <- psl_oof_loghaz_dt(
#     oof_buffer = oof_buffers[[k]],
#     z_covariates = z_covariates
#   )
#
#   dev_k <- poisson_deviance_by_folder_cols(
#     log_hazard_cols = loghaz_cols,
#     tij = as.numeric(dt_by_cause_cv[[k]][["tij"]]),
#     delta = as.integer(dt_by_cause_cv[[k]][["deltaij"]]),
#     folder = as.integer(dt_by_cause_cv[[k]][["folder"]]),
#     nfold = nfold
#   )
#
#   dev_sum <- dev_sum + dev_k
# }
#
#   dev_mean <- colMeans(2.0 * dev_sum)
#
#   dt_learners <- data.table::data.table(
#     learner = learners_labels,
#     deviance = dev_mean
#   )

  ## ------------------------------------------------------------
  ## Step 4: after CV, remove only learners that are ALL-NA
  ## within at least one cause. Partial fold failures are allowed.
  ## ------------------------------------------------------------




  # failed_by_risk <- lapply(oof_buffers, function(buf) {
  #   z_covariates[colSums(!is.na(buf$Z)) == 0L]
  # })
  #
  # failed_cv_learners <- Reduce(union, failed_by_risk)
  #
  # if (length(failed_cv_learners) > 0L) {
  #   keep_z <- setdiff(z_covariates, failed_cv_learners)
  #   keep_ix <- match(keep_z, z_covariates)
  #
  #   oof_buffers <- lapply(oof_buffers, function(buf) {
  #     buf$Z <- buf$Z[, keep_ix, drop = FALSE]
  #     buf
  #   })
  #
  #   learners <- learners[keep_ix]
  #   learners_labels <- learners_labels[keep_ix]
  #   z_covariates <- keep_z
  #   full_train_list <- lapply(full_train_list, function(fits_k) fits_k[keep_ix])
  #
  #   dt_learners <- dt_learners[learner %in% learners_labels]
  # }

  if (length(z_covariates) == 0L) {
    stop("All learners failed: every cross-validated prediction column was entirely NA in at least one competing risk.")
  }

  ## If post-CV pruning leaves one learner, return direct fit
  if (length(learners) == 1L) {
    message("Only one usable base learner remains after cross-validation screening. No ensemble constructed.")

    one_learner_out <- vector("list", n_crisks)

    for (cause_ix in seq_len(n_crisks)) {
      one_learner_out[[cause_ix]] <- list(
        learners_fit = full_train_list[[cause_ix]][[1L]],
        meta_learner_fit = NULL
      )
    }

    out <- list(
      learners = learners,
      metalearner = NULL,
      superlearner = one_learner_out,
      cross_validation_deviance = dt_learners,
      data_info = list(
        id = id,
        status = status,
        event_time = event_time,
        nodes = sort(unique(grid_nodes)),
        nfold = nfold,
        maximum_followup = maximum_followup,
        n_crisks = n_crisks,
        learners_labels = learners_labels,
        variable_transformation = variable_transformation
      )
    )

    class(out) <- "poisson_superlearner"
    return(out)
  }

  ## ------------------------------------------------------------
  ## Step 5: meta-learning, without duplicating the meta-learner
  ## ------------------------------------------------------------


  # meta_learner <- Learner_glmnet(
  #   covariates = z_covariates,
  #   cross_validation = FALSE,
  #   intercept = FALSE,
  #   add_nodes = FALSE,
  #   penalise_nodes = TRUE,
  #   lambda = 0
  # )
  #
  # meta_learner_fits <- mapply(
  #   function(dt_k, oof_k) {
  #     fit_meta_learner(
  #       dt = dt_k,
  #       dt_z = oof_k,
  #       meta_learner = meta_learner,
  #       z_covariates = z_covariates
  #     )
  #   },
  #   dt_by_cause,
  #   oof_buffers,
  #   SIMPLIFY = FALSE
  # )
  #
  # superlearner_out <- lapply(seq_len(n_crisks), function(k) {
  #   list(
  #     learners_fit = full_train_list[[k]],
  #     meta_learner_fit = meta_learner_fits[[k]]
  #   )
  # })

  out <- list(
    learners = learners,
    metalearner = meta_learner,
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
      learners_labels = learners_labels,
      variable_transformation = variable_transformation
    )
  )

  class(out) <- "poisson_superlearner"
  out
}
