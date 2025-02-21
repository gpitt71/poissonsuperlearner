learner_nn <- function(data,covariates){


  return(0)



}

learner_xgboost <- function(data,covariates){


  return(0)



}



#' \code{glmnet} learner
#'
#' Learner that uses the \code{glmnet} algorithm.
#'
#' @param data \code{data.table}, input data to be pre-processed.
#' @param id \code{character}, identifier column.
#' @param cross_validation \code{logical}, whether the learner hyper-parameters need to be learnt for each data split.
#' @param learner_label \code{character} and optional, user-definer label for the learner.
#' @param ... any additional argument to pass to \code{glmnet} or \code{cv.glmnet}.
#'
#' @return \code{data.table} containing the pre-processed data.
#'
#'
#' @export
learner_glmnet <- function(covariates,
                           treatment=NULL,
                           cross_validation=FALSE,
                           ...){

  # browser()
  xs <- paste(covariates, collapse="+")

  if(!is.null(treatment)){

    xs <- paste(xs,"+",treatment)
  }

  formula_string <- paste("deltaij ~", xs,"-1+node+offset(log(tij))", sep="")

  if(cross_validation){

    learner <- function(data,
                        formula=formula_string,
                        ...){

      tmp <- datapp_glmnet(data, formula)

      fit <- cv.glmnet(tmp$x,
                       tmp$y,
                       offset = tmp$offset,
                       family="poisson", ...)

      return (fit)
      }

  }else{

    learner <- function(data,
                        formula=formula_string,
                        ...){

      tmp <- datapp_glmnet(data, formula)

      fit <- glmnet(tmp$x,
                       tmp$y,
                       offset = tmp$offset,
                       family="poisson", ...)

      return (fit)
    }

  }


  return(learner)



}
