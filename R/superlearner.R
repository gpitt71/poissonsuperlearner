#' Fit a Poisson Super Learner ensemble
#'
#' Fits an ensemble of cause-specific hazard learners in long-format Poisson
#' regression and combines them through a meta-learner.
#'
#' @param data `data.frame`. Input subject-level data.
#' @param id `character(1)`. Name of the subject identifier column. If missing,
#'   an `id` column is created automatically.
#' @param status `character(1)`. Name of the event-status column (`0` for censoring,
#'   positive integers for causes).
#' @param event_time `character(1)`. Name of the event/censoring time column.
#' @param learners `list`. List of initialized learner reference-class objects
#'   (e.g. [Learner_glmnet()], [Learner_hal()], [Learner_gam()]).
#' @param number_of_nodes `numeric(1)` or `NULL`. Number of quantile-based
#'   time nodes. Ignored when `nodes` is supplied.
#' @param nodes `numeric` or `NULL`. Explicit time grid for piecewise-constant
#'   hazards.
#' @param variable_transformation `list`/`character`/`formula` or `NULL`.
#'   Transformation(s) evaluated in the long-format data (e.g.
#'   `"x2 ~ I(x^2)"`).
#' @param nfold `numeric(1)`. Number of folds for cross-validation stacking.
#' @param ... Additional arguments currently ignored.
#'
#' @return An object of class `poisson_superlearner` (list) with components
#'   `learners`, `metalearner`, `superlearner`, and `data_info`.
#'
#' @examples
#' data <- simulateStenoT1(200, competing_risks = TRUE)
#' learners <- list(
#'   glm = Learner_glmnet$new(covariates = c("age", "value_LDL"), add_nodes = TRUE),
#'   gam = Learner_gam$new(covariates = c("age", "value_LDL"))
#' )
#' fit <- Superlearner(
#'   data = data, id = "id", status = "status_cvd", event_time = "time_cvd",
#'   learners = learners, number_of_nodes = 10, nfold = 3
#' )
#'
#' @export
Superlearner <- function(data,
                         id = "id", #
                         status = "status", #
                         event_time = NULL, #
                         learners,
                         number_of_nodes = NULL, #
                         nodes = NULL,
                         variable_transformation = NULL,
                         nfold = 3, #
                         ...) {


    if (!(id %in% names(data))) {
        data[["id"]] <- 1:NROW(data)
        id <- "id"
    }
  # give names to the learners
  if(is.null(names(learners))){

    names(learners) <- paste0("learner_",1:length(learners))
  }

  learners_labels = names(learners)

  # Data pre-processing ----

  # save some relevant values

  maximum_followup = max(data[[event_time]])

  n <- length(unique(data[[id]]))#nrow(data)

  min_depth = length(learners)


  if (!(0 %in% data[[status]])) {
    warning(
      paste0(
        "There is no value of ",
        status,
        " equal to zero. We will consider the data uncensored."
      )
    )
    n_crisks <- length(unique(data[[status]]))
    uncensored_01 <-TRUE

  } else{
    n_crisks <- length(unique(data[[status]])) - 1
    uncensored_01 <-FALSE
  }


    # Data pre-processing ----

    #  Handle nodes
    ##Either the nodes are given or we take all of the realised times
    if (!is.null(number_of_nodes)) {

      grid_nodes = quantile(data[[event_time]], probs = seq(0, 1, length.out = as.integer(number_of_nodes)+1), type = 1, names = FALSE)



    } else{
      if (is.null(nodes)) {
        grid_nodes <- sort(unique(data[[event_time]]))

      } else{
        grid_nodes <- nodes

      }
    }

    # Add zero if missing
    if (!(0 %in% grid_nodes)) {
      grid_nodes <- c(0, grid_nodes)

    }



    grid_nodes <- grid_nodes[grid_nodes <= max(data[[event_time]])]

    # Actual data pp
    dt <- data_pre_processing(
      data = data,
      id = id,
      status = status,
      nodes = grid_nodes,
      event_time = event_time,
      uncensored_01=uncensored_01
    )






  lhs_string = NULL


  ## Transform the variables if needed ----
  if (!is.null(variable_transformation)) {apply_transformations(dt,variable_transformation)}

  ## Splitting in folds ----

  id_fold <- sample(1:nfold,
                      n,
                      replace = TRUE,
                      prob = rep(1 / nfold, nfold))

  dt_id <- data.table(folder = id_fold, id = unique(data[[id]]))

  dt <- merge(dt, dt_id, by = "id", all.x = T)


  dt_z <- vector("list", n_crisks)

  z_covariates <- paste0("Z", 1:length(learners))

  # if only one learner is present, we simply perform a CV ----
  if (length(learners) == 1) {

    message("Only one base learner supplied. Fitting the learner directly; no ensemble constructed.")

    # The learner on the full dataset ----
    training_data <- split(dt, by = "k")
    learner_fit <- mapply(function(x){

      out <- learners[[1]]$private_fit(x)
      return(out)
    },
    training_data,
    SIMPLIFY = FALSE
      )

    one_learner_out <- list()

    for(causes in 1:n_crisks){


      one_learner_out[[causes]] <- list(
        model = NULL,
        learners_fit=learner_fit[[causes]],
        meta_learner_fit = NULL
      )

    }

    out <- list(
      learners = learners,
      metalearner = NULL, # it does not exict in this scenario
      superlearner = one_learner_out,
      data_info = list(
        id = id,
        status = status,
        event_time = event_time,
        nodes = sort(unique(as.numeric(levels(dt$node)))),
        nfold = nfold,
        maximum_followup = maximum_followup,
        n_crisks=n_crisks,
        variable_transformation = variable_transformation
      )
    )

    class(out) <- "poisson_superlearner"

    return(out)

  }


  # Train your models in the training set ----
  for (ix in 1:nfold) {


    # Training data for each competing risk ----
    tmp_train <- dt[folder != ix, ]
    training_data <- split(tmp_train, by = "k")
    # Validation data for each competing risk ----
    tmp_val <- dt[folder == ix, ]
    validation_data <- split(tmp_val, by = "k")

    # we find the pseudo observations for each fold ----
    pseudo_observations <- mapply(
      function(training_data,
               validation_data,
               competing_risk,
               learners,
               z_covariates,
               ix)
        create_pseudo_observations(training_data, validation_data, competing_risk, learners, z_covariates, ix),
      training_data = training_data,
      validation_data = validation_data,
      competing_risk=as.list(1:n_crisks),
      MoreArgs = list(
        learners = learners,
        ix = ix,
        z_covariates = z_covariates
      ),
      SIMPLIFY = FALSE
    )



    dt_z <- mapply(function(x, y)
      rbind(x, y), dt_z, pseudo_observations, SIMPLIFY = FALSE)

  }


  L <- length(z_covariates)

  dev_sum <- matrix(0.0, nrow = nfold, ncol = L)

  for (k in seq_len(n_crisks)) {


    loghaz_cols <- (dt_z[[k]][, ..z_covariates])

    dev_k <- poisson_deviance_by_folder_cols(
      log_hazard_cols = loghaz_cols,
      tij    = as.numeric(dt_z[[k]][["tij"]]),
      delta  = as.integer(dt_z[[k]][["deltaij"]]),
      folder = as.integer(dt_z[[k]][["folder"]]),
      nfold  = nfold
    )
    dev_sum <- dev_sum + dev_k
  }

  # 2 * sum_i(...) per fold, then mean across folds
  dev_mean <- colMeans(2.0 * dev_sum)

  dt_learners <- data.table::data.table(
    learner  = paste0("learner_", seq_len(L)),
    deviance = dev_mean
  )



  z_covariates_list <- select_covariate_path(dt_learners, z_covariates, min_depth = min_depth)
  ##

  ## Meta learning

  meta_learner <- Learner_glmnet(
    covariates = z_covariates,
    cross_validation = FALSE,
    intercept = FALSE,
    add_nodes = FALSE,
    penalise_nodes = TRUE,
    lambda=0
  )

  meta_learner_fits <- mapply(
    function(dt,
             dt_z,
             meta_learner,
             learners,
             z_covariates)
      fit_meta_learner(dt, dt_z, meta_learner, learners, z_covariates),
    split(dt, by = "k"),
    dt_z,
    MoreArgs = list(
      meta_learner = meta_learner,
      learners = learners,
      z_covariates = z_covariates
    ),
    SIMPLIFY = FALSE


  )


  dt_learners[,covariate:=NULL]

  dt_learners[,learner:=learners_labels]

  out <- list(
    learners = learners,
    metalearner = meta_learner,
    superlearner = meta_learner_fits,
    cross_validation_deviance = dt_learners,
    data_info = list(
      id = id,
      status = status,
      event_time=event_time,
      nodes = sort(unique(as.numeric(levels(dt$node)))),
      nfold = nfold,
      maximum_followup = maximum_followup,
      n_crisks=n_crisks,
      learners_labels=learners_labels,
      variable_transformation = variable_transformation
    )
  )

  class(out) <- "poisson_superlearner"

  return(out)



}
