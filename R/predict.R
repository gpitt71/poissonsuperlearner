#' Individual Data Pre-Processing
#'
#' This function predicts the survivor function or the hazard using the superlearner.
#'
#' @return predictions
#'
#' @export
predict.poisson_superlearner <- function(object,
                                         newdata,
                                         times,
                                         type = "survival",
                                         cause = 1,
                                         ...) {
  setDT(newdata)

  tmp <- copy(newdata)

  # here we disregard the event_time column if present in the newdata
  tmp <- tmp[, setdiff(names(tmp), object$data_info$event_time), with = FALSE]


  cond_zero <- 0 %in% times

  cond_times_larger_than_max <- times > object$data_info$maximum_followup


  if(all(cond_times_larger_than_max)){

    warning(paste0("All the entries in the input times are larger than the maximum follow-up: ",
                   as.character(object$data_info$maximum_followup)))
    d <- NULL

  }else{

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

  # browser()

  # if (is.null(data_pp[[object$data_info$id]])) {
    # data_pp[[object$data_info$id]] <- 1:nrow(data_pp)
  # }

  # no problem writing over id
  data_pp[[object$data_info$id]] <- 1:nrow(data_pp)


  # if (is.null(data_pp[[object$data_info$status]])) {
    data_pp[[object$data_info$status]] <- 0
  # }

  data_pp <- data_pre_processing(
    data_pp,
    id = object$data_info$id,
    status = object$data_info$status,
    event_time = object$data_info$event_time,
    nodes = object$data_info$nodes
  )



  # browser()
  if(object$data_info$matrix_transformation){

    # browser()

    columns_of_interest <- unlist(lapply(object$learners, function(x){return(unique(c(x$covariates,x$treatment)))}))

    columns_of_interest <- unique(columns_of_interest[(complete.cases(columns_of_interest))])

    # Take the variable that we transform

    lhs_vars <- trimws(unlist(strsplit(strsplit(object$data_info$variable_transformation, "~")[[1]][1], "\\+")))
    lhs_string <- paste(lhs_vars, collapse = ", ")

    # Take the transformation
    rhs_vars <- trimws(unlist(strsplit(strsplit(object$data_info$variable_transformation, "~")[[1]][2], "\\+")))
    rhs_string <- paste(rhs_vars, collapse = ", ")


    eval(parse(text = paste0(
      "
               data_pp[,c('", lhs_string

      , "'):=list(", rhs_string
      , ")]
               "
    )))


    # data_pp <- data_pp[,.(tij = sum(tij),
    #             deltaij=sum(deltaij)), by = c(unique(c(columns_of_interest,lhs_string)),"node","k")]
    #
    #
    # data_pp[,c("id"):=1:nrow(data_pp)]



  }


  # Set covariates for metalearner
  z_covariates <- paste0("Z", 1:length(object$learners))


  # Predict on the validation set your pseudo-observations ----
  # browser()

  data_pp[,deltatime:=tij][,tij:=1]

  if(length(object$learners)==1){


    # browser()

    # dt_pred <- object$learners[[1]]$predictor(
    #   model=object$superlearner$learners_fit[[1]],
    #   newdata = data_pp
    #
    # )

    dt_pred <- mapply(
      function(crisk_cause, model,superl_fit, newdata)
        model$predictor(model = superl_fit$learners_fit, newdata = data_pp),
      as.list(1:object$data_info$n_crisks),
     object$superlearner,
      MoreArgs = list(newdata = data_pp, model=object$learners[[1]]),
     SIMPLIFY = F
    )


  }else{





  learners_predictions <- mapply(
    function(f, model, newdata)
      f$predictor(model = model, newdata = data_pp),
    object$learners,
    object$superlearner[[cause]]$learners_fit,
    MoreArgs = list(newdata = data_pp)
  )

  # browser()

  pseudo_observations_data <- matrix(apply(as.matrix(learners_predictions, nrow=nrow(newdata), ncol=length(z_covariates)), MARGIN = 2, log),
                                     nrow=nrow(data_pp),
                                     ncol=length(z_covariates))



  # Name the columns

  colnames(pseudo_observations_data) <- z_covariates

  setDT(as.data.frame.matrix(pseudo_observations_data))

  dt_pred <- object$superlearner[[cause]]$model$predictor(object$superlearner[[cause]]$meta_learner_fit,
                                                          newdata =
                                                            cbind(pseudo_observations_data, data_pp))


  }

  # browser()

  data_pp[,paste0("pwch_",1:object$data_info$n_crisks):=dt_pred]

  pwch_cols <- paste0("pwch_",1:object$data_info$n_crisks)

  mapply(function(pwch, name) {
    data_pp[, (paste0("cumulative_hazard_", name)) := cumsum(get(pwch) * deltatime), by = id]
  }, pwch_cols, gsub("pwch_", "", pwch_cols))


  hazard_terms <- paste0("cumulative_hazard_", 1:object$data_info$n_crisks)
  sum_expr <- paste(hazard_terms, collapse = " + ")
  survival_function_string <- paste0("data_pp[, survival_function := exp(-(", sum_expr, "))]")

  eval(parse(text = survival_function_string))

  # data_pp[,survival_function_shift := shift(survival_function,fill=1),by=id]

  # absolute_risk_string <- paste0(
  #   "data_pp[, absolute_risk := cumsum(survival_function_shift * pwch_", cause, " * deltatime), by = id]"
  # )

  absolute_risk_string <- paste0(
    "data_pp[, absolute_risk := cumsum(survival_function * pwch_", cause, " * deltatime), by = id]"
  )

  eval(parse(text = absolute_risk_string))

  # data_pp <- copy(data_pp)
  # data_pp[,survival_function:=pmin(exp(-cumsum(pwch*deltatime)),1), by=.(id)]

  data_pp <- data_pp[, .SD[.N], by = id][,times:=as.numeric(as.character(node))+deltatime]


  if(cond_zero){

    data_pp[time==0,
            c(pwch_cols,
              'survival_function'):=list(rep(0,length(pwch_cols)),1)]



  }

  columns_ss <- unique(c(colnames(newdata),object$data_info$event_time,pwch_cols,"survival_function","absolute_risk"))

  d <- data_pp[,..columns_ss]

  }



  if(any(cond_times_larger_than_max)){

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
    d2[,c(pwch_cols,
            'survival_function'):=list(rep(NA,length(pwc_cols)),NA)]


    if (object$data_info$id %in% colnames(d)) {
      d2[[object$data_info$id]] <- (nrow(data_pp)+1):(nrow(data_pp)+nrow(d2))
    }



    d <- rbind(d,d2)


  }


  return(d)




}
