#' Poisson Super Learner
#'
#' This class is used for estimating the Poisson SuperLearner.
#'
#' @param data \code{data.frame}, input data to be pre-processed.
#' @param id \code{character}, identifier column.
#' @param status \code{character}, status column.
#' @param event_time \code{character}, time-to-event column in competing risks and survival applications.
#' @param learners \code{list}, list of learners to include in the ensemble. If only one learner is provided, the learner is applied to the full data.
#' @param variable_transformation \code{list}, variable transformation(s) to apply to the data. Each element of the list is a character \code{"new_variable ~ some_function(variable_in_the_data)"}
#' @param nfold \code{numeric}, number of V-folds to construct the ensemble.
#' @param number_of_nodes \code{numeric}, number of time points sampled from the observed time to construct the nodes. Alternative to \code{nodes}, if both NULL we take all the observed time points as nodes.
#' @param nodes \code{numeric}, time grid to construct the piece-wise constant model. Alternative to \code{number_of_nodes}, if both NULL we take all the observed time points as nodes.
#'
#' @return A \code{poisson_superlearner} object contains the following output
#' \itemize{
#' \item{\code{learners}: \code{list} containing the learners.}
#' \item{\code{metalearner}: \code{numeric} Training negative log likelihood.}
#' \item{\code{superlearner}:  \code{numeric} Validation  negative log likelihood. Not available for COX.}
#' \item{\code{meta_learner_cross_validation}:  \code{data.table} containing the average Poisson Deviance evaluated on the V-Folds for the learners and the meta-learner.}
#' \item{\code{data_info}: \code{list} containing the following meta-data:}
#'    \itemize{
#'    \item{\code{id}: \code{character}, name of the covariate in the data that contains the id variable.}
#'    \item{\code{status}: \code{character}, name of the covariate in the data that contains the status variable.}
#'    \item{\code{event_time}: \code{character}, for competing risks or survival data it contains the name of the covariate containing the event time.}
#'    \item{\code{nfold}: \code{integer} denoting the number of V-Folds.}
#'    \item{\code{maximum_followup}: \code{numeric} denoting the maximum follow-up time observed.}
#'    \item{\code{n_crisks}, \code{integer} denoting the number of competing risks.}
#'    \item{\code{variable_transformation}: \code{list} or \code{array}, containing the data variable transformations implemented. }
#'    \item{\code{interval_data_type}: \code{logical}, \code{TRUE} for interval data. }
#'    }
#' \item{\code{hazard_model}: \code{string} chosen hazard model (COX, NN or XGB)}
#' \item{\code{IndividualDataPP}: starting \code{IndividualDataPP} object.}
#' }
#' @export
Superlearner <- function(data,
                         id = "id", #
                         status = "status", #
                         event_time = NULL, #
                         learners,
                         number_of_nodes = NULL, #
                         nodes = NULL,
                         min_depth=NULL,
                         meta_learner_algorithms = c("glm","glmnet"),
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

    maximum_followup = max(data[[event_time]])

  # Data pre-processing ----

  # save some relevant values
  n <- length(unique(data[[id]]))#nrow(data)




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


  data_by_competing_risk <- split(dt, by = "k")

  # Compute oos deviance scores for the learners ----
  dt_learners <- rbindlist(
    lapply(seq_along(dt_z), function(i) {
      DT <- copy(dt_z[[i]])
      # detect all Z* columns (flexible: Z1, Z_2, Zabc3, etc.)
      zcols <- grep("^Z", names(DT), value = TRUE)
      melt(
        DT,
        id.vars = c("id", "folder", "node"),
        measure.vars = zcols,
        variable.name = "learner",
        value.name = "pwch"
      )[, which := paste0("pwch_", i)]
    }),
    use.names = TRUE, fill = TRUE
  )

  # transform back to exponential
  dt_learners[,pwch:=exp(pwch)]

  dt_learners[, learner_idx := fifelse(
    grepl("\\d+", learner),
    as.integer(gsub("\\D+", "", learner)),
    match(learner, unique(learner))  # stable ordering for non-numeric suffixes
  )]
  dt_learners[, learner := paste0("learner_", learner_idx)][, learner_idx := NULL]

  dt_learners <- dcast(
    dt_learners,
    id + folder + node + learner ~ which,
    value.var = "pwch"
    # , fun.aggregate = mean  # uncomment if duplicates exist within a cell
  )

  dt_learners[,node:=factor(node,levels = sort(as.numeric(levels(node))),ordered=TRUE)]

  setorder(dt_learners, id, node, learner)
  dt_learners <- dt_learners[complete.cases(dt_learners),]

  # Actual deviance computation

  pwch_cols <- paste0("pwch_",1:n_crisks)

  # save sum of pwch

  sum_of_hazards <- paste(pwch_cols, collapse = " + ")

  pwch_dot_string <- paste0("dt_learners[, pwch_dot :=",sum_of_hazards,"]")

  eval(parse(text = pwch_dot_string))

  dt_learners<- merge(dt_learners, dt[k==1,.(id,node,tij)], by = c("id", "node"))

  # compute cumulative hazard

  mapply(function(pwch, name) {
    dt_learners[, (paste0("cumulative_hazard_", name)) := cumsum(get(pwch) * tij), by = .(id,learner)]
  }, pwch_cols, gsub("pwch_", "", pwch_cols))


  # compute survival function

  hazard_terms <- paste0("cumulative_hazard_", 1:n_crisks)
  sum_expr <- paste(hazard_terms, collapse = " + ")
  survival_function_string <- paste0("dt_learners[, survival_function := exp(-(", sum_expr, "))]")

  eval(parse(text = survival_function_string))

  mapply(function(pwch, name) {
    dt_learners[, (paste0("hazard_times_time_", name)) := (get(pwch) * tij)]
  }, pwch_cols, gsub("pwch_", "", pwch_cols))


  # get the deltas
  delta_list =mapply(function(risk,data){data[k==risk,][,c("id","node",paste0("delta_",risk)):=list(id,node,deltaij)][,.SD,.SDcols=c("id","node",paste0("delta_",risk))]},as.list(1:n_crisks),MoreArgs=list(data=dt),SIMPLIFY = F)

  delta_list <- merge_deltas(delta_list)

  delta_list[,node:=factor(node,levels = sort(as.numeric(levels(node))),ordered=TRUE)]

  dt_learners<- merge(dt_learners,
                      delta_list,
                      by=c("id","node"),
                      all.x = T)


  #first term

  mapply(function(name) {
    dt_learners[,paste0("cr_contribution_",name):=get(paste0("delta_", name))*log(get(paste0("delta_", name))/get(paste0("hazard_times_time_", name)))]
    dt_learners[,paste0("cr_contribution_",name):=fifelse(is.nan(get(paste0("cr_contribution_",name))),0,get(paste0("cr_contribution_",name)))]
    dt_learners[,paste0("cr_contribution_",name):=get(paste0("cr_contribution_",name))-(get(paste0("delta_", name))-get(paste0("hazard_times_time_", name)))]

  }, as.list(1:n_crisks))


  # second term

  crc_terms <- paste0("cr_contribution_", 1:n_crisks)
  sum_expr <- paste(crc_terms, collapse = " + ")
  crc_string <- paste0("dt_learners[, cr_contribution_tot := ",sum_expr,"]")

  eval(parse(text = crc_string))


  dt_learners<-dt_learners[,.(deviance_i=sum(cr_contribution_tot)),by=.(id,learner)]

  dt_learners <- merge(dt_learners,dt_id,by="id")

  setkey(dt_learners,NULL)

  dt_learners<-dt_learners[, .(deviance_v = 2*sum(deviance_i, na.rm = TRUE)), by = .(learner,folder)]

  dt_learners<-dt_learners[,.(deviance=mean(deviance_v)),by =learner]


  z_covariates_list <- select_covariate_path(dt_learners, z_covariates, min_depth = min_depth)


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
    data_by_competing_risk,
    dt_z,
    MoreArgs = list(
      meta_learner = meta_learner,
      learners = learners,
      z_covariates = z_covariates
    ),
    SIMPLIFY = FALSE


  )



  setnames(dt_learners,"learner","model")

  dt_learners[,covariate:=NULL]

  out <- list(
    learners = learners,
    metalearner = meta_learner,
    superlearner = meta_learner_fits,
    # meta_learner_cross_validation=rbind(dt_cv_out,
    #                                     dt_learners),
    data_info = list(
      id = id,
      status = status,
      event_time=event_time,
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
