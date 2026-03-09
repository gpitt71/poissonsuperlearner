#' Fit a Poisson Super Learner ensemble
#'
#' Fits an ensemble of **cause-specific** piecewise-constant hazard models using a
#' long-format Poisson representation and combines them through a **meta-learner**
#' (stacking). Cross-validation is used to obtain out-of-sample base-learner
#' predictions for meta-learning and to compute cross-validated deviances for the
#' base learners.
#'
#' Internally, the function:
#' \enumerate{
#'   \item Builds a time grid (`nodes`) and converts the subject-level data to a long
#'   Poisson format (one row per subject-interval-cause).
#'   \item Splits subjects into `nfold` folds.
#'   \item For each fold and each cause, fits each base learner on the training folds
#'   and predicts on the held-out fold to create out-of-sample pseudo-observations
#'   (columns `Z1`, `Z2`, ...).
#'   \item Fits a (cause-specific) meta-learner on the stacked pseudo-observations.
#'   \item Refits each base learner on the full long data (per cause) and stores all
#'   fitted objects for prediction.
#' }
#'
#' @param data `data.frame`. Subject-level input data (one row per subject).
#' @param id `character(1)`. Name of the subject identifier column. If missing,
#'   an `id` column is created automatically.
#' @param status `character(1)`. Name of the event-status column.
#'   Must be coded with `0` = censoring and `1,2,...,K` for event types (causes).
#'   If there is no `0` in `status`, the data are treated as uncensored.
#' @param event_time `character(1)`. Name of the event/censoring time column.
#' @param learners `list`. List of initialized learner reference-class objects
#'   (e.g. [Learner_glmnet()], [Learner_hal()], [Learner_gam()]). If unnamed, learners
#'   are named `"learner_1"`, `"learner_2"`, ...
#'   Each learner must implement `$private_fit(dt_long)` and `$private_predictor(model, newdata)`.
#' @param number_of_nodes `numeric(1)` or `NULL`. If not `NULL`, constructs a
#'   quantile-based node grid with `number_of_nodes + 1` cut points (including
#'   endpoints), then adds `0` if missing. Ignored when `nodes` is supplied.
#' @param nodes `numeric` or `NULL`. Explicit time-node grid (cut points). If
#'   supplied, `number_of_nodes` is ignored. `0` is added if missing. Nodes beyond
#'   `max(event_time)` are dropped.
#' @param variable_transformation `list`/`character`/`formula` or `NULL`.
#'   Optional transformations applied to the internally created long Poisson data
#'   before fitting (via `apply_transformations()`).
#' @param nfold `numeric(1)`. Number of folds for cross-validation stacking.
#' @param ... Additional arguments currently ignored.
#'
#' @return An object of class `poisson_superlearner`, i.e. a named `list` with:
#'
#' \describe{
#' \item{learners}{The list of base learner objects as provided in `learners`
#'   (possibly auto-named). These store each learner specification (covariates,
#'   tuning parameters, etc.).}
#'
#' \item{metalearner}{The meta-learner object used to combine base learners.
#'   Currently this is a `Learner_glmnet(...)` created internally with
#'   `covariates = c("Z1","Z2",...)`, `intercept = FALSE`, and `lambda = 0`.
#'   \strong{Learner-dependent details:} only the meta-learner \emph{object} is stored
#'   here; the \emph{fitted} meta-learner for each cause is in `superlearner[[k]]$meta_learner_fit`.
#'
#'   If only one base learner is supplied, `metalearner` is `NULL` and no ensemble is
#'   constructed (see below).}
#'
#' \item{superlearner}{A `list` of length `data_info$n_crisks` (one entry per cause).
#'   Entry `superlearner[[k]]` (cause `k`) is itself a list with:
#'   \describe{
#'     \item{model}{The meta-learner object (same class as `metalearner`).}
#'     \item{learners_fit}{A `list` of fitted base-learner model objects (one per
#'       element of `learners`), refit on the \emph{full} long data for cause `k`.
#'       Each element is \strong{learner-dependent} (e.g. `"glmnet"`/`"fishnet"`,
#'       `"gam"`, etc.) and is the object passed to that learner’s
#'       `$private_predictor()` method at prediction time.}
#'     \item{meta_learner_fit}{The fitted meta-learner object for cause `k`, fit on the
#'       stacked long data containing the out-of-sample base-learner predictions
#'       (`Z1`, `Z2`, ...) and the original long data columns. This object is
#'       \strong{learner-dependent} (for the default meta-learner it is a `"glmnet"`
#'       fit).}
#'   }
#'
#'   \strong{Single-learner special case:} if `length(learners) == 1`, then
#'   `superlearner[[k]]$model` and `superlearner[[k]]$meta_learner_fit` are `NULL`,
#'   and `superlearner[[k]]$learners_fit` contains the fitted base-learner model for
#'   cause `k` (no stacking performed).}
#'
#' \item{cross_validation_deviance}{A `data.table` with columns:
#'   \describe{
#'     \item{learner}{Learner labels (taken from `names(learners)`).}
#'     \item{deviance}{Mean cross-validated Poisson deviance for each base learner,
#'       aggregated across causes and averaged over folds.}
#'   }
#'   This is only present when `length(learners) > 1`.}
#'
#' \item{data_info}{A `list` of bookkeeping information needed for prediction and
#'   interpretation:
#'   \describe{
#'     \item{id}{Identifier column name used.}
#'     \item{status}{Status column name used.}
#'     \item{event_time}{Event/censoring time column name used.}
#'     \item{nodes}{Numeric vector of node cut points used for the piecewise grid
#'       (includes `0` and is sorted). These are the interval boundaries used in the
#'       long Poisson representation.}
#'     \item{nfold}{Number of folds used for stacking.}
#'     \item{maximum_followup}{`max(data[[event_time]])`.}
#'     \item{n_crisks}{Number of event types (causes) detected.
#'       If censoring is present (`0` in `status`), then `n_crisks = #unique(status) - 1`;
#'       otherwise `n_crisks = #unique(status)`.}
#'     \item{learners_labels}{Character vector of learner labels (names of the `learners` list).}
#'     \item{variable_transformation}{The transformation specification passed in
#'       `variable_transformation` (or `NULL`).}
#'   }}
#' }
#'
#' @examples
#' data <- simulateStenoT1(200, competing_risks = TRUE)
#' learners <- list(
#'   glm = Learner_glmnet(covariates = c("age", "value_LDL"), lambda = 0, cross_validation = FALSE),
#'   gam = Learner_gam(covariates = c("age", "value_LDL"))
#' )
#' fit <- Superlearner(
#'   data = data, id = "id", status = "status_cvd", event_time = "time_cvd",
#'   learners = learners, number_of_nodes = 10, nfold = 3
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

  if (is.null(names(learners))) {
    names(learners) <- paste0("learner_", seq_along(learners))
  }

  learners_labels <- names(learners)

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
    grid_nodes <- nodes
  }

  if (!(0 %in% grid_nodes)) {
    grid_nodes <- c(0, grid_nodes)
  }

  grid_nodes <- grid_nodes[grid_nodes <= max(data[[event_time]])]

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

  z_covariates <- paste0("Z", seq_along(learners))

  ## ------------------------------------------------------------
  ## Step 1: fit learners once on the full data, one list per cause
  ## ------------------------------------------------------------
  full_train_list <- lapply(dt_by_cause, function(dt_k) {
    lapply(learners, function(f) f$private_fit(dt_k))
  })

  ## Remove learners that already fail on the full data
  failed_by_full_fit <- lapply(full_train_list, function(fits_k) {
    z_covariates[vapply(fits_k, is_failed_fit, logical(1))]
  })

  failed_learners <- Reduce(union, failed_by_full_fit)

  if (length(failed_learners) > 0L) {
    keep_z <- setdiff(z_covariates, failed_learners)
    keep_ix <- match(keep_z, z_covariates)

    learners <- learners[keep_ix]
    learners_labels <- learners_labels[keep_ix]
    z_covariates <- keep_z

    full_train_list <- lapply(full_train_list, function(fits_k) fits_k[keep_ix])
  }

  if (length(z_covariates) == 0L) {
    stop("All learners failed on the full data before cross-validation.")
  }

  ## If pruning leaves one learner, return direct fit
  if (length(learners) == 1L) {
    message("Only one usable base learner remains. Fitting the learner directly; no ensemble constructed.")

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
      data_info = list(
        id = id,
        status = status,
        event_time = event_time,
        nodes = sort(unique(as.numeric(levels(dt$node)))),
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
  ## Step 2: V-fold CV only on retained learners
  ## ------------------------------------------------------------
  dt_z <- vector("list", n_crisks)

  for (ix in seq_len(nfold)) {
    tmp_train <- dt[folder != ix, ]
    training_data <- split(tmp_train, by = "k")

    tmp_val <- dt[folder == ix, ]
    validation_data <- split(tmp_val, by = "k")

    pseudo_observations <- mapply(
      function(training_data,
               validation_data,
               competing_risk,
               learners,
               z_covariates,
               ix) {
        create_pseudo_observations(
          training_data = training_data,
          validation_data = validation_data,
          competing_risk = competing_risk,
          learners = learners,
          z_covariates = z_covariates,
          ix = ix
        )
      },
      training_data = training_data,
      validation_data = validation_data,
      competing_risk = as.list(seq_len(n_crisks)),
      MoreArgs = list(
        learners = learners,
        ix = ix,
        z_covariates = z_covariates
      ),
      SIMPLIFY = FALSE
    )

    dt_z <- mapply(function(x, y) rbind(x, y), dt_z, pseudo_observations, SIMPLIFY = FALSE)
  }

  ## ------------------------------------------------------------
  ## Step 3: CV deviance on retained learners
  ## ------------------------------------------------------------

  L <- length(z_covariates)
  dev_sum <- matrix(0.0, nrow = nfold, ncol = L)

  for (k in seq_len(n_crisks)) {
    loghaz_cols <- dt_z[[k]][, ..z_covariates]

    dev_k <- poisson_deviance_by_folder_cols(
      log_hazard_cols = loghaz_cols,
      tij = as.numeric(dt_z[[k]][["tij"]]),
      delta = as.integer(dt_z[[k]][["deltaij"]]),
      folder = as.integer(dt_z[[k]][["folder"]]),
      nfold = nfold
    )

    dev_sum <- dev_sum + dev_k
  }

  dev_mean <- colMeans(2.0 * dev_sum)

  dt_learners <- data.table::data.table(
    learner = learners_labels,
    deviance = dev_mean
  )

  ## ------------------------------------------------------------
  ## Step 4: after CV, remove only learners that are ALL-NA
  ## within at least one cause. Partial fold failures are allowed.
  ## ------------------------------------------------------------
  failed_by_risk <- lapply(dt_z, function(DT) {
    z_covariates[vapply(DT[, ..z_covariates], function(x) all(is.na(x)), logical(1))]
  })

  failed_cv_learners <- Reduce(union, failed_by_risk)

  if (length(failed_cv_learners) > 0L) {
    keep_z <- setdiff(z_covariates, failed_cv_learners)
    keep_ix <- match(keep_z, z_covariates)

    dt_z <- lapply(dt_z, function(DT) {
      DT[, c(setdiff(names(DT), z_covariates), keep_z), with = FALSE]
    })

    learners <- learners[keep_ix]
    learners_labels <- learners_labels[keep_ix]
    z_covariates <- keep_z
    full_train_list <- lapply(full_train_list, function(fits_k) fits_k[keep_ix])

    dt_learners <- dt_learners[learner %in% learners_labels]
  }

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
        nodes = sort(unique(as.numeric(levels(dt$node)))),
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
  meta_learner <- Learner_glmnet(
    covariates = z_covariates,
    cross_validation = FALSE,
    intercept = FALSE,
    add_nodes = FALSE,
    penalise_nodes = TRUE,
    lambda = 0
  )

  meta_learner_fits <- mapply(
    function(dt_k, dt_z_k) {
      fit_meta_learner(
        dt = dt_k,
        dt_z = dt_z_k,
        meta_learner = meta_learner,
        z_covariates = z_covariates
      )
    },
    dt_by_cause,
    dt_z,
    SIMPLIFY = FALSE
  )

  superlearner_out <- lapply(seq_len(n_crisks), function(k) {
    list(
      learners_fit = full_train_list[[k]],
      meta_learner_fit = meta_learner_fits[[k]]
    )
  })

  out <- list(
    learners = learners,
    metalearner = meta_learner,
    superlearner = superlearner_out,
    cross_validation_deviance = dt_learners,
    data_info = list(
      id = id,
      status = status,
      event_time = event_time,
      nodes = sort(unique(as.numeric(levels(dt$node)))),
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
