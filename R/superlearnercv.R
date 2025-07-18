#' K-fold Cross-Validation Poisson Super Learner
#'
#' This function implements the Poisson Super Learner.
#'
#' @param data \code{data.frame}, input data to be pre-processed.
#' @param id \code{character}, identifier column.
#' @param status \code{character}, status column.
#' @param nodes \code{numeric}, time grid to construct the piece-wise constant model.
#'
#' @return superlearner
#'
#' @export
Superlearner.cv <- function(data,
                         id = "id",
                         stratified_k_fold=FALSE,
                         start_time = NULL,
                         end_time = NULL,
                         status = "status",
                         event_time = NULL,
                         learners,
                         n_cross_validation_folds=5,
                         number_of_nodes = NULL,
                         nodes = NULL,
                         meta_learner_algorithm = "glmnet",
                         add_nodes_metalearner=TRUE,
                         add_intercept_metalearner=TRUE,
                         matrix_transformation=FALSE,
                         nested_cross_validation_meta_learner=TRUE,
                         penalise_nodes_metalearner=TRUE,
                         variable_transformation=NULL,
                         nfold = 3,
                         ...) {



  # browser()

  # save some relevant values
  n <- length(unique(data[[id]]))#nrow(data)
  n_crisks <- length(unique(data[[status]])) - 1


  if (stratified_k_fold) {

    setDT(data)

    dt_id <- eval(parse(text = paste0("data[,last(",status,"),by = ",
                                      id,"]")))

    eval(parse(text = paste0("setnames(dt_id,'V1','",status,"')")))

    dt_id<-stratified_sampling(dt_id,id,status,n_cross_validation_folds)

    dt_id[order(id)]


  } else{
    id_fold <- sample(1:n_cross_validation_folds,
                      n,
                      replace = TRUE,
                      prob = rep(1 / n_cross_validation_folds, n_cross_validation_folds))

    dt_id <- data.table(folder = id_fold, id = unique(data[[id]]))
  }





  dt <- merge(data,dt_id,by="id")


  likelihood_train <- likelihood_val <- NULL

  for(n_cv_f in 1:n_cross_validation_folds){

    tmp_train <- dt[folder != n_cv_f, ][,folder:=NULL]
    # Validation data for each competing risk ----
    tmp_val <- dt[folder == n_cv_f, ][,folder:=NULL]

    sl_model <- Superlearner(tmp_train,
                             id = id,
                             stratified_k_fold=stratified_k_fold,
                                            start_time = start_time,
                                            end_time = end_time,
                                            status = status,
                                            event_time = event_time,
                                            learners,
                                            number_of_nodes=number_of_nodes,
                                            nodes = nodes,
                                            meta_learner_algorithm = meta_learner_algorithm,
                                            add_nodes_metalearner=add_nodes_metalearner,
                                            add_intercept_metalearner=add_intercept_metalearner,
                                            matrix_transformation=matrix_transformation,
                                            penalise_nodes_metalearner=penalise_nodes_metalearner,
                                            variable_transformation=variable_transformation,
                                            nfold = nfold,
                             nested_cross_validation_meta_learner=nested_cross_validation_meta_learner,
                             ...)




    # browser()
    # evaluate in parallell all the cause-specific super-learners ----


    train_ev <- lapply(X=1:n_crisks,
                       FUN=cv_subject_specific_hazard,
                       object=sl_model,
                            newdata=tmp_train)

    validation_ev <- lapply(X=1:n_crisks,
                       FUN=cv_subject_specific_hazard,
                       object=sl_model,
                       newdata=tmp_val)

      # cv_subject_specific_hazard()


    combined_train_ev <- rbindlist(train_ev)
    combined_val_ev <- rbindlist(validation_ev)


    likelihood_train <- c(likelihood_train,
                          combined_train_ev[,likelihood_i:=sum(deltaij*log(pwch)-log(pwch)),by=id][,sum(likelihood_i)])

    likelihood_val <- c(likelihood_val,
                        combined_val_ev[,likelihood_i:=sum(deltaij*log(pwch)-log(pwch)),by=id][,sum(likelihood_i)])


  }



  # browser()

  out_cv <- data.table(

    stratified_k_fold=stratified_k_fold,
    n_nodes = min(number_of_nodes,
                  length(nodes)),
    meta_learner_algorithm = meta_learner_algorithm,
    add_nodes_metalearner=add_nodes_metalearner,
    add_intercept_metalearner=add_intercept_metalearner,
    matrix_transformation=matrix_transformation,
    penalise_nodes_metalearner=penalise_nodes_metalearner,
    variable_transformation=variable_transformation,
    nfold = nfold,
    train_likelihood=mean(likelihood_train,
                          na.rm=TRUE)/1000,
    validation_likelihood=mean(likelihood_val,
                               na.rm=TRUE)/1000
  )

  message("Likelihood values are divided by 1000.")

  return(out_cv)

}
