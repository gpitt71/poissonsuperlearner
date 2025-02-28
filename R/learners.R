#' Script containing our learners
#'
#' This script contains the learners available in the package.
#'
#' @import data.table

learner_nn <- function(data,covariates){


  return(0)



}

learner_xgboost <- function(data,covariates){


  return(0)



}

#' \code{glm} learner class
#'
#' @export Learner_glm
#' @exportClass Learner_glm
Learner_glm <- setRefClass(
  "Learner_glm",
  fields = list(
    covariates = "character",
    treatment = "character",
    cross_validation = "logical",
    competing_risks ="logical",
    formula ="character",
    learner="function"
  ),
  methods = list(
    initialize = function(covariates = NULL,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          competing_risks = NA ,
                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$cross_validation <- cross_validation

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula(covariates = .self$covariates,
                                      treatment = .self$treatment,
                                      competing_risks =TRUE)


      if (.self$cross_validation) {
        .self$learner = NULL

      } else{
        .self$learner = glm
      }
    },

    fit = function(data, ...) {

      survival_01 <- length(unique(data[['k']])) < 2


      # practical correction in case it is a survival problem
      if(survival_01){

        .self$formula <- create_formula(covariates = .self$covariates,
                                        treatment = .self$treatment,
                                        competing_risks =FALSE)

      }

      out<- .self$learner(.self$formula,
                          data=data,
                          # offset = tmp$offset,
                          family = "poisson",
                          ...)

      return(out)

    },

    predictor = function(model, newdata,...) {

      # tmp <- datapp_glmnet(newdata, .self$formula)

      out <- predict(model,
                     ...,
                     newdata=newdata,
                     type = "response")

      return(out)


    }



  )
)



#' \code{glmnet} learner class
#'
#' @export Learner_glmnet
#' @exportClass Learner_glmnet
Learner_glmnet <- setRefClass(
  "Learner_glmnet",
  fields = list(
    covariates = "character",
    treatment = "character",
    cross_validation = "logical",
    # competing_risks="logical",
    formula ="character",
    learner="function"
  ),
  methods = list(

    initialize = function(covariates = NULL,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          # competing_risks = NA,
                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$cross_validation <- cross_validation

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula(covariates = .self$covariates,
                                      treatment = .self$treatment,
                                      competing_risks =TRUE)


      if (.self$cross_validation) {
        .self$learner = cv.glmnet

      } else{
        .self$learner = glmnet
      }
    },

    fit = function(data, formula, ...) {

      survival_01 <- length(unique(data[['k']])) < 2


      # practical correction in case it is a survival problem
      if(survival_01){

        .self$formula <- create_formula(covariates = .self$covariates,
                                        treatment = .self$treatment,
                                        competing_risks =FALSE)

      }

      # one layer of extra data preprocessing
      tmp <- datapp_glmnet(data, .self$formula)

      out<- .self$learner(tmp$x,
                           tmp$y,
                           offset = tmp$offset,
                           family = "poisson",
                           ...)

      return(out)

    },

    predictor = function(model, newdata,...) {

      tmp <- datapp_glmnet(newdata, .self$formula)

      out <- predict(model,
                     ...,
                   newx=tmp$x,
                   newoffset = tmp$offset,
                   type = "response")

      return(out)


    }



  )
)


# learner_glmnet <- function(covariates,
#                            treatment=NULL,
#                            cross_validation=FALSE,
#                            ...){
#
#   # browser()
#   xs <- paste(covariates, collapse="+")
#
#   if(!is.null(treatment)){
#
#     xs <- paste(xs,"+",treatment)
#   }
#
#   formula_string <- paste("deltaij ~", xs,"-1+node+offset(log(tij))", sep="")
#
#   if(cross_validation){
#
#     learner <- function(data,
#                         formula=formula_string,
#                         ...){
#
#       tmp <- datapp_glmnet(data, formula)
#
#       fit <- cv.glmnet(tmp$x,
#                        tmp$y,
#                        offset = tmp$offset,
#                        family="poisson", ...)
#
#       return (fit)
#     }
#
#   }else{
#
#     learner <- function(data,
#                         formula=formula_string,
#                         ...){
#
#       tmp <- datapp_glmnet(data, formula)
#
#       fit <- glmnet(tmp$x,
#                     tmp$y,
#                     offset = tmp$offset,
#                     family="poisson", ...)
#
#       return (fit)
#     }
#
#   }
#
#
#   out <- list(learner=learner,
#               predictor)
#
#
#   class(out) <- "learner_glmnet"
#
#   return(out)
#
#
#
# }














