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

  if (type != "loss") {
    tmp <- copy(newdata)
    # here we disregard the event_time column if present in the newdata
    tmp <- tmp[, setdiff(names(tmp), object$data_info$event_time), with = FALSE]


    cond_zero <- 0 %in% times

    cond_times_larger_than_max <- times > object$data_info$maximum_followup


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

      tmp[, dummy := 1]
      vec_dt[, dummy := 1]

      # Merge on dummy to create Cartesian product
      data_pp <- merge(tmp, vec_dt, by = "dummy", allow.cartesian = TRUE)[, dummy := NULL]
    }
  } else{

    data_pp <- copy(newdata)

    status_per_id <- newdata[, .(status = max(get(object$data_info$status))), by = c(object$data_info$id)]


  }

  # browser()

  # if (is.null(data_pp[[object$data_info$id]])) {
    # data_pp[[object$data_info$id]] <- 1:nrow(data_pp)
  # }

  # no problem writing over id
  if (is.null(data_pp[[object$data_info$id]])) {
    data_pp[[object$data_info$id]] <- 1:nrow(data_pp)
  }

  if (is.null(data_pp[[object$data_info$status]])) {
    data_pp[[object$data_info$status]] <- 0
  }

  data_pp <- data_pre_processing(
    data_pp,
    id = object$data_info$id,
    status = object$data_info$status,
    event_time = object$data_info$event_time,
    nodes = object$data_info$nodes
  )



  # browser()
  if(!is.null(object$data_info$variable_transformation)){

    # browser()

    # columns_of_interest <- unlist(lapply(object$learners, function(x){return(unique(c(x$covariates,x$treatment)))}))

    # columns_of_interest <- unique(columns_of_interest[(complete.cases(columns_of_interest))])

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


    #

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


    # browser()

    learners_predictions <- mapply(
      function(crisk_cause,superlearner,newdata,learners)
        learners_hat(crisk_cause,superlearner,newdata,learners),
      crisk_cause=as.list(1:object$data_info$n_crisks),
      superlearner=object$superlearner,
      MoreArgs = list(newdata = data_pp,
                      learners=object$learners),
      SIMPLIFY = FALSE
    )


    pseudo_observations_data <- lapply(
      learners_predictions,
      function(mx){
        out <- matrix(apply(
        as.matrix(
          mx, nrow=nrow(newdata), ncol=length(z_covariates)), MARGIN = 2, log),
        nrow=nrow(data_pp),
        ncol=length(z_covariates))
        # out <-as.matrix(mx)
        colnames(out) <- z_covariates # Name the columns
        as.data.table(as.data.frame.matrix(out))}
    )


  dt_pred <- mapply(
    function(crisk_cause,superlearner,pseudo_observations_data){
      superlearner$model$predictor(superlearner$meta_learner_fit,newdata =cbind(pseudo_observations_data, data_pp))},
    as.list(1:object$data_info$n_crisks),
    object$superlearner,
    pseudo_observations_data,
    SIMPLIFY = F
  )


  }

  # browser()

  # save casue-specific pwch

  data_pp[,paste0("pwch_",1:object$data_info$n_crisks):=dt_pred]

  pwch_cols <- paste0("pwch_",1:object$data_info$n_crisks)

  # save sum of pwch

  sum_of_hazards <- paste(pwch_cols, collapse = " + ")

  pwch_dot_string <- paste0("data_pp[, pwch_dot :=",sum_of_hazards,"]")

  eval(parse(text = pwch_dot_string))


  # compute cumulative hazard

  mapply(function(pwch, name) {
    data_pp[, (paste0("cumulative_hazard_", name)) := cumsum(get(pwch) * deltatime), by = id]
  }, pwch_cols, gsub("pwch_", "", pwch_cols))


  # compute survival function

  hazard_terms <- paste0("cumulative_hazard_", 1:object$data_info$n_crisks)
  sum_expr <- paste(hazard_terms, collapse = " + ")
  survival_function_string <- paste0("data_pp[, survival_function := exp(-(", sum_expr, "))]")

  eval(parse(text = survival_function_string))

  # shift survival function

  data_pp[,survival_function_shift := shift(survival_function,fill=1),by=id]

  absolute_risk_string <- paste0(
    "data_pp[, absolute_risk := cumsum(survival_function_shift * pwch_", cause, "/pwch_dot * (1-exp(-pwch_dot*deltatime))), by = id]"
  )


  eval(parse(text = absolute_risk_string))


  # this is essentially the likelihood computation
  if(type=="loss"){

  lkh_dt <- copy(data_pp)

  lkh_dt <- merge(lkh_dt, status_per_id, by = "id", all.x = TRUE)

  cols_delta <- paste0("delta_", 1:object$data_info$n_crisks)

  lkh_dt<-lkh_dt[, (cols_delta) := lapply(seq_along(cols_delta), function(j) {
    # j corresponds to delta_j
    if (unique(status) == j) {
      # only the column matching status gets last row = 1
      out <- rep(0L, .N)
      out[.N] <- 1L
      out
    } else {
      # all other columns are zeros
      rep(0L, .N)
    }
  }), by = id]

  mapply(function(pwch, name) {
    lkh_dt[, (paste0("hazard_times_time_", name)) := (get(pwch) * deltatime)]
  }, pwch_cols, gsub("pwch_", "", pwch_cols))

  lkh_dt[,paste0("cr_contribution_",cause):=deltaij*log(get(paste0("hazard_times_time_", cause)))-get(paste0("hazard_times_time_", cause))]


  other_causes <- gsub("pwch_", "", pwch_cols)
  other_causes <- other_causes[other_causes != cause]


  mapply(function(name) {
    lkh_dt[,paste0("cr_contribution_",name):=get(paste0("delta_", name))*log(get(paste0("hazard_times_time_", name)))-get(paste0("hazard_times_time_", name))]
  }, other_causes)



  crc_terms <- paste0("cr_contribution_", 1:object$data_info$n_crisks)
  sum_expr <- paste(crc_terms, collapse = " + ")
  crc_string <- paste0("lkh_dt[, cr_contribution_tot := ",sum_expr,"]")

  eval(parse(text = crc_string))


  lkh_dt<-lkh_dt[,.(log_likelihood_i=sum(cr_contribution_tot)),by=id]

  return(lkh_dt)

}


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
