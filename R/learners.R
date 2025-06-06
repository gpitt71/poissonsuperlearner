#' Script containing our learners
#'
#' This script contains the learners available in the package.
#'
#' @import data.table

Learner_nn <- function(data,covariates){


  return(0)



}


Learner_xgboost <- function(data,covariates){


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
                          intercept = TRUE,
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

    # fit = function(data, validation_data=NULL, ...) {

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




    #   out<- .self$learner(.self$formula,
    #                       data=rbind(data,
    #                                  validation_data),
    #                       family = "poisson",
    #                       subset = rep(TRUE,dim(data)[1]))
    #
    #   return(out)
    #
    # },


    fit = function(data, ...) {

      # out<- .self$learner(.self$formula,
      #                     data=data,
      #                     family = "poisson", ...)

      # ss <- data[['train_01']]
      out <- glm.fit(x=sparse.model.matrix(formula(.self$formula),
                                           data),#[ss,]
                     y=data[['deltaij']], #[ss]
                     offset = log(data[['tij']]), #[ss]
                     family = poisson())

      return(out)

    },


    private_fit = function(data, ...) {

      # browser()
      # ss <- data[['train_01']]
      out <- glm.fit(x=sparse.model.matrix(formula(.self$formula),
                                           data),#[ss,],
                     y=data[['deltaij']],#[ss],
                     offset = log(data[['tij']]),#[ss]),
                     family = poisson())

      # out<- .self$learner(.self$formula,
      #                     data=data[train_01==TRUE,],
      #                     family = "poisson", ...)

      return(out)


    },

    predictor = function(model, newdata,...) {

      newdata<-newdata[complete.cases(newdata),]

      newx=sparse.model.matrix(formula(.self$formula),
                               newdata)

      newoffset = log(newdata[['tij']])

      # browser()

      out <-exp(newx %*% model$coefficients + newoffset)

      # out <- predict(model,
      #                ...,
      #                newdata=newdata,
      #                type = "response")

      return(as.numeric(out))
      # return(out)


    },

    private_predictor = function(model, newdata, ...) {



      newdata<-newdata[complete.cases(newdata),]

      # ss <- !newdata[['train_01']]


      # out <- predict(model,
      #                ...,
      #                newdata=newdata[train_01==FALSE,],
      #                type = "response")

      # return(out)



      newx=sparse.model.matrix(formula(.self$formula),
                               newdata)#[ss,]

      newoffset = log(newdata[['tij']])#[ss])

      # browser()

      out <-exp(newx %*% model$coefficients + newoffset)


      return(as.numeric(out))






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
    recycle_information = "logical",
    intercept="logical",
    formula ="character",
    learner="function",
    add_nodes="logical",
    penalise_nodes= "logical",
    fit_arguments = "list",
    covariates_attributes_matrix= "list"
  ),
  methods = list(

    initialize = function(covariates = NULL,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          intercept=FALSE,
                          add_nodes = TRUE,
                          penalise_nodes=FALSE,
                          recycle_information =FALSE,

                          ...) {
      .self$covariates <- covariates

      .self$treatment <- treatment

      .self$cross_validation <- cross_validation

      .self$intercept <- intercept

      .self$add_nodes <- add_nodes

      .self$penalise_nodes <- penalise_nodes

      .self$recycle_information <- recycle_information

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

      .self$covariates_attributes_matrix <- list(...)

      .self$fit_arguments[['family']] <- "poisson"

      .self$fit_arguments[['intercept']] <- .self$intercept




    },


    datapp = function(train_data=NULL,
                      validation_data=NULL) {





      if(!is.null(train_data)){


        train.mf  <- model.frame(as.formula(.self$formula),
                                 rbind(train_data,
                                       validation_data),
                                 drop.unused.levels = FALSE)

        x  <- model.matrix(attr(train.mf, "terms"), data = train_data)

        y  <- train_data[['deltaij']]

        offset <- log(train_data[['tij']])

        .self$recycle_information <- TRUE

        .self$covariates_attributes_matrix[['train.mf']] <- train.mf



      }else{

        x  <- model.matrix(attr(.self$covariates_attributes_matrix[['train.mf']], "terms"), data = validation_data)
        y  <- validation_data[['deltaij']]
        offset <- log(validation_data[['tij']])


      }


      out <- list(x = x, y = y, offset = offset)


      return(out)
    },

    fit = function(data, ...) {

      data<-data[complete.cases(data),]

      x = sparse.model.matrix(formula(.self$formula),
                              data)

      .self$fit_arguments[['y']] <- as.numeric(data[['deltaij']])

      .self$fit_arguments[['offset']] <-  log(data[['tij']])


      if(!.self$penalise_nodes){

        .self$fit_arguments[['penalty.factor']] <- 1- grepl("node",colnames(x))

      }

      .self$fit_arguments[['x']] <- x

      out <- do.call(.self$learner,
                     .self$fit_arguments)

      return(out)

    },

    private_fit = function(data, ...) {


      # browser()

      # ss <- data[['train_01']]
      x = sparse.model.matrix(formula(.self$formula),
                              data)#[ss,]

      if(!.self$penalise_nodes){

        .self$fit_arguments[['penalty.factor']] <- 1- grepl("node",colnames(x))

      }




      .self$fit_arguments[['x']] <- x
      .self$fit_arguments[['y']] <- data[['deltaij']]#[ss]
      .self$fit_arguments[['offset']] <- log(data[['tij']])#[ss])


      out <- do.call(.self$learner,
                     .self$fit_arguments)


      return(out)

    },

    predictor = function(model, newdata, ...) {


      # if(.self$recycle_information){
      #
      #   tmp <- .self$datapp(validation_data = newdata)
      #
      #   }else{
      #
      #     tmp <- datapp_glmnet(newdata, .self$formula)
      #
      #   }

      # browser()
      out <- predict(model,
                     ...,
                   newx=sparse.model.matrix(formula(.self$formula),
                                            newdata),
                   newoffset = log(newdata[['tij']]),
                   type = "response")

      return(out)


    },

    private_predictor = function(model, newdata, ...) {

      out <- predict(model,
                     ...,
                     newx=sparse.model.matrix(formula(.self$formula),
                                              newdata),#[ss,],
                     newoffset = log(newdata[['tij']]),#[ss]),
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














