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
    intercept ="logical",
    formula ="character",
    learner="function",
    add_nodes="logical"
  ),
  methods = list(
    initialize = function(covariates = NULL,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          intercept = FALSE,
                          add_nodes= TRUE,
                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$cross_validation <- cross_validation

      .self$add_nodes <- add_nodes

      .self$intercept<-intercept

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula(covariates = .self$covariates,
                                      treatment = .self$treatment,
                                      intercept = .self$intercept,
                                      add_nodes=.self$add_nodes)


      if (.self$cross_validation) {
        .self$learner = glm
        warning("There is no cross-validation procedure for a glm model. This learner is a simple glm.")

      } else{
        .self$learner = glm
      }
    },

    fit = function(data, ...) {

      # survival_01 <- length(unique(data[['k']])) < 2


      # practical correction in case it is a survival problem
      # if(survival_01){
      #
      #   .self$formula <- create_formula(covariates = .self$covariates,
      #                                   treatment = .self$treatment,
      #                                   competing_risks =FALSE,
      #                                   add_nodes=.self$add_nodes)
      #
      # }

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
    intercept="logical",
    formula ="character",
    learner="function",
    add_nodes="logical",
    fit_arguments = "list"
  ),
  methods = list(

    initialize = function(covariates = NULL,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          intercept=FALSE,
                          add_nodes = TRUE,
                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$cross_validation <- cross_validation

      .self$intercept <- intercept

      .self$add_nodes <- add_nodes

      # create formula for competing risks. It is correct in the fit method if survival.
      .self$formula <- create_formula(covariates = .self$covariates,
                                      treatment = .self$treatment,
                                      intercept =FALSE, #in the glmnet case, intercept is handled separately.
                                      add_nodes=.self$add_nodes)


      if (.self$cross_validation) {
        .self$learner = cv.glmnet

      } else{
        .self$learner = glmnet
      }

      .self$fit_arguments <- list(...)

      .self$fit_arguments[['family']] <- "poisson"

      .self$fit_arguments[['intercept']] <- .self$intercept


    },

    fit = function(data, formula) {

      # survival_01 <- length(unique(data[['k']])) < 2


      # practical correction in case it is a survival problem
      # if(survival_01){
      #
      #   .self$formula <- create_formula(covariates = .self$covariates,
      #                                   treatment = .self$treatment,
      #                                   competing_risks =FALSE,
      #                                   add_nodes=.self$add_nodes)
      #
      # }

      # one layer of extra data preprocessing
      tmp <- datapp_glmnet(data, .self$formula)

      .self$fit_arguments[['x']] <- tmp$x
      .self$fit_arguments[['y']] <- tmp$y
      .self$fit_arguments[['offset']] <- tmp$offset


      out <- do.call(.self$learner,
                     .self$fit_arguments)


      # out<- .self$learner(,
      #                      tmp$y,
      #                      offset = ,
      #                      family = "poisson",
      #                     intercept = FALSE,
      #                      ...)

      return(out)

    },

    predictor = function(model, newdata, ...) {

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














