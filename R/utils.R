# Data preprocessing ----

create_response_variable_c_risks <- function(nodes, time_to_event, delta, event_type){


  p_holder <- ifelse(delta == event_type, 1, 0)

  l <- sum(nodes < time_to_event)

  out <- c(rep(0, max(0, l - 1)),
           p_holder)

  return(out)
}

create_offset_variable <- function(nodes, delta, time_to_event){

  # if(time_to_event==0){
  #
  #   return(cbind(0,0))
  #
  # }else{
  tmp <- c(nodes[nodes < time_to_event],
           first(nodes[nodes >= time_to_event]))



  if (all(nodes < time_to_event)) {
    tmp <- c(tmp, time_to_event)
  } else{
    tmp[length(tmp)] <- time_to_event
  }


  tij <- diff(c(tmp))

  grid_nodes <- c(nodes[nodes < time_to_event])

  return(cbind(grid_nodes,tij))


  # }
}

data_pre_processing <- function(data,
                                id,
                                status,
                                event_time,
                                nodes=NULL,
                                predictions=FALSE,
                                uncensored_01=FALSE
){


  setDT(data)

  # Handle competing risks ----
  ## for each of the competing risks (CR) we need to create a table
  n_crisks <- pmax(length(unique(data[[status]])) - 1+uncensored_01,1)
  ## the CR tables are stuck on top of each other to allow for possible interactions
  dt_fit <- do.call(rbind, replicate(n_crisks, data, simplify = FALSE))
  ## we create an artificial k index. Table specific.
  dt_fit <- dt_fit[, k := rep(1:n_crisks, each = dim(data)[1])]

    # Data Transformation ----
    tmp <- c(id, "k")
    dt_fit <- dt_fit[ , {
        tte <- .SD[[1L]]   # time-to-event vector
        del <- .SD[[2L]]   # status/delta vector
        off <- create_offset_variable(nodes, time_to_event = tte)
        .(node    = off[, 1L],
          tij     = off[, 2L],
          deltaij = create_response_variable_c_risks(nodes,time_to_event = tte,delta = del,event_type = k)
          )
    }, by = c(id, "k"),
    .SDcols = c(event_time, status)]

    ## Retrieve covariates

  dt_fit <- merge(dt_fit, data, by = id, all.x = TRUE)

  setnames(dt_fit, c(id),c("id"))

  if(predictions){
    maxn <-last(nodes)
    lvls <- as.character(sort(unique(nodes)))


  }else{

    maxn <- max(dt_fit$node)
    lvls <- as.character(sort(unique(dt_fit$node)))
  }



  dt_fit[,c("node",
            "k"):=list(factor(node, levels=lvls),
                       as.factor(k))]

  dt_fit[,node:=relevel(node,ref=as.character(maxn))]

  # dt_fit[,node:=relevel(node,ref=as.character(last(nodes)))]

  return(dt_fit)

}



## Matrix transformation ----

apply_transformations <- function(dt, variable_transformation) {
  stopifnot(data.table::is.data.table(dt))

  if (is.null(variable_transformation)) return(invisible(dt))

  # Normalize to a list of items (strings or formulas)
  items <- variable_transformation
  if (inherits(items, "formula")) items <- list(items)
  else if (is.character(items))    items <- as.list(items)
  else if (is.list(items))         items <- items
  else stop("Unsupported 'variable_transformation' type.")

  # Collect all (lhs name, rhs expression) pairs
  lhs_all  <- character()
  rhs_all  <- vector("list", 0L)

  for (it in items) {
    if (inherits(it, "formula")) {
      # formula: lhs ~ rhs
      lhs_raw <- paste(deparse(it[[2L]]), collapse = "")
      rhs_raw <- paste(deparse(it[[3L]]), collapse = "")
    } else if (is.character(it) && length(it) == 1L) {
      parts <- strsplit(it, "~", fixed = TRUE)[[1L]]
      if (length(parts) != 2L) stop("Each transformation must contain a single '~'.")
      lhs_raw <- parts[1L]
      rhs_raw <- parts[2L]
    } else {
      stop("List elements must be formulas or single strings of the form 'lhs ~ rhs'.")
    }

    lhs_vec <- trimws(strsplit(lhs_raw, "\\+")[[1L]])
    rhs_vec <- trimws(strsplit(rhs_raw, "\\+")[[1L]])

    if (length(lhs_vec) != length(rhs_vec)) {
      stop(
        sprintf("LHS and RHS have different lengths in '%s ~ %s' (%d vs %d).",
                lhs_raw, rhs_raw, length(lhs_vec), length(rhs_vec))
      )
    }

    # Parse each RHS into an expression
    rhs_exprs <- lapply(rhs_vec, function(s) parse(text = s)[[1L]])

    lhs_all <- c(lhs_all, lhs_vec)
    rhs_all <- c(rhs_all, rhs_exprs)
  }

  # Evaluate and assign each transformation.
  # Using one-by-one assignment keeps evaluation within data.table's environment
  # so column names are visible and functions come from the calling env.
  for (i in seq_along(lhs_all)) {
    dt[, (lhs_all[i]) := eval(rhs_all[[i]])]
  }

  invisible(dt)
}

sl_cut <- function(x, breaks, include.lowest = TRUE, right = TRUE) {

  labels <- paste0(head(breaks, -1), "-", tail(breaks, -1) - if (right) 1 else 0)
  # Return labeled factor
  cut(x, breaks = breaks, labels = labels, include.lowest = include.lowest, right = right)
}


# Other utils ----
create_formula <- function(covariates=NA_character_,
                           competing_risks=FALSE,
                           intercept=FALSE,
                           add_nodes=TRUE){


 xs<- NULL

  if (!any(is.na( covariates))) {
  xs <- paste(covariates, collapse = "+")
  }

  # if (!is.na( treatment)) {
  #   xs <- paste(xs, "+", treatment)
  # }

  # if (competing_risks) {
  #   xs <- paste(xs, "+ k")
  # }

  if (add_nodes) {
    xs <- paste(xs, "+ node")
  }

  if(!intercept){
    xs <- paste(xs, "-1")
  }

  out <- paste("deltaij ~", xs, "+offset(log(tij))", sep =
                 "")

  return(out)




}


create_formula_glmnet <- function(covariates=NA_character_,
                           add_nodes=TRUE){


  xs<- NULL

  if (!any(is.na( covariates))) {
    xs <- paste(covariates, collapse = "+")
  }

  # if (!is.na( treatment)) {
  #   xs <- paste(xs, "+", treatment)
  # }

  if (add_nodes) {
    xs <- paste(xs, "+ node")
  }

  out <- paste("deltaij ~", xs, "+offset(log(tij))", sep =
                 "")

  return(out)




}

create_formula_gam <- function(covariates=NA_character_,
                           competing_risks=FALSE,
                           intercept=FALSE,
                           add_nodes=TRUE){



  if (!any(is.na( covariates))) {
    xs <- paste(covariates, collapse = "+")
  }

  # if (!is.na( treatment)) {
  #   xs <- paste(xs, "+", treatment)
  # }

  # if (competing_risks) {
  #   xs <- paste(xs, "+ k")
  # }

  if (add_nodes) {
    xs <- paste(xs, "+ node")
  }

  # if(!intercept){
  #   xs <- paste(xs, "-1")
  # }

  out <- paste("deltaij ~", xs, sep =
                 "")

  return(out)




}

create_formula_hal <- function(covariates=NA_character_,
                           competing_risks=FALSE,
                           intercept=FALSE,
                           add_nodes=TRUE){



  if (!any(is.na( covariates))) {
    xs <- paste(covariates, collapse = "*")
  }

  # if (!is.na( treatment)) {
  #   xs <- paste(xs, "*", treatment)
  # }

  # if (competing_risks) {
  #   xs <- paste(xs, "+ k")
  # }

  if (add_nodes) {
    xs <- paste(c(xs, "node"), collapse = "*")
  }

  if(!intercept){
    xs <- paste(xs, "-1")
  }

  out <- paste("deltaij ~", xs, "+offset(log(tij))", sep =
                 "")

  return(out)




}



create_pseudo_observations <- function(training_data,
                                       validation_data,
                                       competing_risk,
                                       learners,
                                       z_covariates,
                                       ix){
  "
  This function creates the pseudo-observations.

  "

  train_list <- lapply(learners, function(f) f$private_fit(training_data))


  # Predict on the validation set your pseudo-observations ----
  val_list <- mapply(
    function(f, model, newdata)
      f$private_predictor(model = model, newdata = newdata),
    learners,
    train_list,
    MoreArgs = list(newdata = validation_data)
  )



  # val_list<- as.matrix(val_list)


  val_list<- apply(as.matrix(val_list),
                 MARGIN = 2,
                 log)


  # Name the columns

  colnames(val_list) <- z_covariates


  dt_z <-  data.table(val_list)[, c("id", "folder","node",paste0("deltaij"), "tij") := validation_data[,.(id,folder,node,deltaij,tij)]] #

  dt_z[,competing_risk:= competing_risk]

  return(dt_z)



}


select_covariate_path <- function(dt, z_covariates, min_depth) {
  # add covariates column
  dt[, covariate := z_covariates]

  # order by deviance
  ordered_cov <- dt[order(deviance), covariate]

  # build list from min_depth up to full length
  lapply(min_depth:length(ordered_cov), function(k) ordered_cov[1:k])
}



merge_deltas <- function(lst) {
  if (length(lst) == 0L) return(data.table(id = integer(), node = integer()))

  # Keep only id, node, and delta_* columns from each table
  cleaned <- lapply(lst, function(DT) {
    DT <- as.data.table(copy(DT))
    dcols <- grep("^delta_", names(DT), value = TRUE)
    if (length(dcols) == 0L) stop("Each element must contain at least one delta_* column.")
    DT <- DT[, c("id", "node", dcols), with = FALSE]
    setkey(DT, id, node)
    DT
  })

  # Full outer merge across the list on id+node
  out <- Reduce(function(x, y) merge(x, y, by = c("id", "node"), all = TRUE), cleaned)

  # Order columns: id, node, then delta_* in natural order (by numeric suffix if present)
  deltas <- grep("^delta_", names(out), value = TRUE)
  ord <- order(suppressWarnings(as.integer(sub(".*?(\\d+)$", "\\1", deltas))), deltas)
  setcolorder(out, c("id", "node", deltas[ord]))

  out[]
}


## Meta learning ----

meta_learners_candidates <- function(meta_learner_algorithms,
                                     z_covariates){

  out <- list()

  ## this should be removed from the production version of the package!!!
  if("glm_ml_1"%in% meta_learner_algorithms){

    one_time_learner=list(glm_ml_1 =  Learner_glmnet(
      covariates = z_covariates,
      cross_validation = FALSE,
      intercept = FALSE,
      add_nodes = FALSE,
      penalise_nodes = TRUE,
      lambda=0
    ))

    out <- c(out,one_time_learner)

  }

  if("glmnet_ml_1"%in% meta_learner_algorithms){

    one_time_learner=list(glmnet_ml_1 =  Learner_glmnet(
      covariates = z_covariates,
      cross_validation = TRUE,
      intercept = FALSE,
      add_nodes = FALSE,
      penalise_nodes = TRUE
    ))

out <- c(out,one_time_learner)

  }

  if("glm_ml_2"%in% meta_learner_algorithms){

    one_time_learner=list(glm_ml_2 =  Learner_glmnet(
      covariates = c(z_covariates, paste0(z_covariates,":node")),
      cross_validation = FALSE,
      intercept = FALSE,
      add_nodes = TRUE,
      penalise_nodes = TRUE,
      lambda=0
    ))

    out <- c(out,one_time_learner)

  }

  if("glmnet_ml_2"%in% meta_learner_algorithms){

    one_time_learner=list(glmnet_ml_2 =  Learner_glmnet(
      covariates = c(z_covariates, paste0(z_covariates,":node")),
      cross_validation = TRUE,
      intercept = FALSE,
      add_nodes = TRUE,
      penalise_nodes = TRUE
    ))

    out <- c(out,one_time_learner)

  }

  if("glm_ml_3"%in% meta_learner_algorithms){

    one_time_learner=list(glm_ml_3 =  Learner_glmnet(
      covariates = z_covariates,
      cross_validation = FALSE,
      intercept = FALSE,
      add_nodes = TRUE,
      penalise_nodes = TRUE,
      lambda=0
    ))

    out <- c(out,one_time_learner)

  }

  if("glmnet_ml_3"%in% meta_learner_algorithms){

    one_time_learner=list(glmnet_ml_3 =  Learner_glmnet(
      covariates =  c(z_covariates, paste0(z_covariates,":node")),
      cross_validation = TRUE,
      intercept = FALSE,
      add_nodes = TRUE,
      penalise_nodes = FALSE
    ))

    out <- c(out,one_time_learner)

  }

  if("glm_ml_4"%in% meta_learner_algorithms){

    one_time_learner=list(glm_ml_4 =  Learner_glmnet(
      covariates =c(z_covariates, paste0(z_covariates,":node")),
      cross_validation = FALSE,
      intercept = FALSE,
      lambda=0,
      add_nodes = FALSE,
      penalise_nodes = TRUE
    ))

    out <- c(out,one_time_learner)

  }

  if("glmnet_ml_4"%in% meta_learner_algorithms){

    one_time_learner=list(glmnet_ml_4 =  Learner_glmnet(
      covariates = z_covariates,
      cross_validation = TRUE,
      intercept = FALSE,
      add_nodes = TRUE,
      penalise_nodes = TRUE
    ))

    out <- c(out,one_time_learner)

  }

  if("glmnet_ml_5"%in% meta_learner_algorithms){

    one_time_learner=list(glmnet_ml_5 =  Learner_glmnet(
      covariates = z_covariates,
      cross_validation = TRUE,
      intercept = FALSE,
      add_nodes = TRUE,
      penalise_nodes = FALSE
    ))

    out <- c(out,one_time_learner)

  }

  if("glmnet_ml_6"%in% meta_learner_algorithms){

    one_time_learner=list(glmnet_ml_6 =  Learner_glmnet(
      covariates =c(z_covariates, paste0(z_covariates,":node")),
      cross_validation = TRUE,
      intercept = FALSE,
      add_nodes = FALSE,
      penalise_nodes = TRUE
    ))

    out <- c(out,one_time_learner)

  }

  #########################################################


  if("glmnet" %in% meta_learner_algorithms){


    glmnet_meta_learners <- list(

      # Z1 + Z2
      glmnet_ml_1 =  Learner_glmnet(
        covariates = z_covariates,
        cross_validation = TRUE,
        intercept = FALSE,
        add_nodes = FALSE,
        penalise_nodes = TRUE
      )#,


      #next paper ----
      ## Z1*Z2*node
      # glmnet_ml_2 =  Learner_glmnet(
      #   covariates = c(z_covariates, paste0(z_covariates,":node")),
      #   cross_validation = TRUE,
      #   intercept = FALSE,
      #   add_nodes = TRUE,
      #   penalise_nodes = TRUE
      # ),
      ## Z1*Z2*node - not penalised
      # glmnet_ml_3 =  Learner_glmnet(
      #   covariates =  c(z_covariates, paste0(z_covariates,":node")),
      #   cross_validation = TRUE,
      #   intercept = FALSE,
      #   add_nodes = TRUE,
      #   penalise_nodes = FALSE
      # ),
      ## Z1+Z2+node
      # glmnet_ml_4 =  Learner_glmnet(
      #   covariates = z_covariates,
      #   cross_validation = TRUE,
      #   intercept = FALSE,
      #   add_nodes = TRUE,
      #   penalise_nodes = TRUE
      # ),
      ## Z1+Z2+node - not penalised
      # glmnet_ml_5 =  Learner_glmnet(
      #   covariates = z_covariates,
      #   cross_validation = TRUE,
      #   intercept = FALSE,
      #   add_nodes = TRUE,
      #   penalise_nodes = FALSE
      # ),
      ## Z1+Z2+Z1:node+Z2:node
      # glmnet_ml_6 =  Learner_glmnet(
      #   covariates =c(z_covariates, paste0(z_covariates,":node")),
      #   cross_validation = TRUE,
      #   intercept = FALSE,
      #   add_nodes = FALSE,
      #   penalise_nodes = TRUE
      # )
      #

      )

    out <- c(out,glmnet_meta_learners)


  }


  if("glm" %in% meta_learner_algorithms){


    glm_meta_learners <- list(

      # Z1 + Z2
      glm_ml_1 =  Learner_glmnet(
        covariates = z_covariates,
        cross_validation = FALSE,
        intercept = FALSE,
        add_nodes = FALSE,
        penalise_nodes = TRUE,
        lambda=0
      )#,
      ## Z1*Z2*node
      # glm_ml_2 =  Learner_glmnet(
      #   covariates = c(z_covariates, paste0(z_covariates,":node")),
      #   cross_validation = FALSE,
      #   intercept = FALSE,
      #   add_nodes = TRUE,
      #   penalise_nodes = TRUE,
      #   lambda=0
      # ),
      ## Z1+Z2+node
      # glm_ml_3 =  Learner_glmnet(
      #   covariates = z_covariates,
      #   cross_validation = FALSE,
      #   intercept = FALSE,
      #   add_nodes = TRUE,
      #   penalise_nodes = TRUE,
      #   lambda=0
      # ),
      ## Z1+Z2+Z1:node+Z2:node
      # glm_ml_4 =  Learner_glmnet(
      #   covariates =c(z_covariates, paste0(z_covariates,":node")),
      #   cross_validation = FALSE,
      #   intercept = FALSE,
      #   lambda=0,
      #   add_nodes = FALSE,
      #   penalise_nodes = TRUE
      # )

      )

    out <- c(out,glm_meta_learners)


  }

  names(out) <- paste(names(out),paste0(z_covariates, collapse = ""),sep = "_")

  return(out)

}


fit_meta_learner <- function(dt,
                             dt_z,
                             meta_learner,
                             learners,
                             z_covariates){

  tmp <- merge(dt_z[, !c("tij", "deltaij"), with = FALSE],dt,by=c("id","folder","node"))

  meta_learner_fit <- meta_learner$private_fit(tmp)

  # learners on the full dataset ----

  full_train_list <- lapply(learners, function(f) f$private_fit(dt))

  out <- list(
    model = meta_learner,
    learners_fit=full_train_list,
    meta_learner_fit = meta_learner_fit
  )



  return(out)


}


# super-learner hyper parameters cross-validation ----


meta_learner_cross_validation <- function(dt,
                             dt_z,
                             cr_ix,
                             nfold,
                             meta_learner){



  # model output ---
  out<-NULL

  dt_z <- merge(dt_z,dt,by=c("id","folder","node"))


  for(v_fold_id in 1:nfold){

  meta_learner_fit <- meta_learner$private_fit(dt_z[folder!=v_fold_id,])

  oos_data <- copy(dt_z[folder == v_fold_id, ])

  # forced tij to be one
  oos_data[,tij:=1]


  fitted_hazard = data.table(
    oos_data[['id']],
    oos_data[['node']],
    oos_data[['deltaij']],
    v_fold_id,
    as.vector(meta_learner$private_predictor(meta_learner_fit, oos_data))
  )

  setnames(fitted_hazard,c("id","node",paste0("delta_",cr_ix),"folder",paste0("pwch_", cr_ix)))

  out <- rbind(out,fitted_hazard)

  }

  return(out)


}



learners_hat <- function(crisk_cause,superlearner,newdata,learners){

  learners_predictions <- mapply(
    function(f, model, newdata)
      f$private_predictor(model = model, newdata = newdata),
    learners,
    superlearner$learners_fit,
    MoreArgs = list(newdata = newdata)
  )

  return(learners_predictions)



}

# Hal helpers ----

mk_main <- function(x, v, K) {
  if (is.numeric(x)) {
    out <- mk_main_numeric_cpp(x, as.integer(K))
    cps <- out$cutpoints
    idxs <- out$idxs

    cps <- as.numeric(cps)
    if (!length(cps)) {
      return(list(idxs = list(), names = character(),
                  prim_meta = list(), var_meta = list(cutpoints = numeric())))
    }

    # nms <- sprintf("I(%s<=%.10g)", v, cps)
    nms <- sprintf("I(%s>=%.10g)", v, cps)

    prim_meta <- lapply(cps, function(cu) list(var = v, kind = "numeric", cutpoint = cu))
    list(
      idxs = idxs,
      names = nms,
      prim_meta = prim_meta,
      var_meta = list(cutpoints = cps)
    )

  } else if (is.factor(x)) {
    lv <- levels(x)
    if (!length(lv)) {
      return(list(idxs = list(), names = character(),
                  prim_meta = list(), var_meta = list(levels = character(), ref = NA_character_)))
    }

    # Drop a reference level (default: first level)
    ref <- lv[1]
    keep <- setdiff(lv, ref)

    # If only one level, nothing to add
    if (!length(keep)) {
      return(list(
        idxs = list(),
        names = character(),
        prim_meta = list(),
        var_meta = list(levels = lv, ref = ref)
      ))
    }

    # Compute idxs for all levels, then drop the reference entry
    out_all <- mk_main_factor_cpp(as.integer(x), length(lv))
    idxs_all <- out_all$idxs

    # Map level -> position (mk_main_factor_cpp is assumed to return in level order)
    pos_keep <- match(keep, lv)

    idxs <- idxs_all[pos_keep]
    nms <- sprintf("I(%s==%s)", v, make.names(keep, unique = TRUE))
    prim_meta <- lapply(keep, function(L) list(var = v, kind = "factor", level = L, ref = ref))

    list(
      idxs = idxs,
      names = nms,
      prim_meta = prim_meta,
      var_meta = list(levels = lv, ref = ref)
    )

  } else {
    stop(sprintf("Unsupported type for %s", v))
  }
}

##
add_cols <- function(state, idxs_list, nm_vec, col_meta_list) {
  out <- add_cols_cpp(idxs_list, p_start = state$p)
  ncol_added <- out$ncol
  if (!ncol_added) return(state)

  keep <- out$keep
  nm_vec <- nm_vec[keep]
  col_meta_list <- col_meta_list[keep]

  state$chunk_used <- state$chunk_used + 1L
  state$I_chunks[[state$chunk_used]] <- out$i
  state$J_chunks[[state$chunk_used]] <- out$j
  state$name_chunks[[state$chunk_used]] <- nm_vec
  state$colmeta_chunks[[state$chunk_used]] <- col_meta_list

  state$p <- state$p + ncol_added
  state
}



