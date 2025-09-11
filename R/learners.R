#' Script containing our learners
#'
#' This script contains the learners available in the package.
#'
#' @import data.table

Learner_nn <- function(data,covariates){


  return(0)



}

#' \code{xgboost} learner class
#'
#' @export Learner_xgboost
#' @exportClass Learner_xgboost
Learner_xgboost <- setRefClass(
  "Learner_xgboost",
  fields = list(
    covariates = "character",
    treatment = "character",
    nrounds = "integer",
    cross_validation = "logical",
    intercept ="logical",
    formula ="character",
    learner="function",
    add_nodes="logical",
    fit_arguments = "list",
    grid_of_hyperparameters = "list"
  ),
  methods = list(
    initialize = function(covariates = NULL,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          intercept = FALSE,
                          add_nodes= TRUE,
                          nrounds = NA_integer_,
                          grid_of_hyperparameters=NULL,
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


      # handle fit arguments


      .self$nrounds<-nrounds

      .self$fit_arguments <- list(...)

      .self$fit_arguments <- list(params=list(objective="count:poisson"))
      .self$fit_arguments[['nrounds']] <-  .self$nrounds

      # handle hyperparamters grid if necessary

      if (.self$cross_validation) {
        .self$grid_of_hyperparameters <- expand.grid(grid_of_hyperparameters)
      } else{

        .self$grid_of_hyperparameters <-grid_of_hyperparameters
        #
        .self$fit_arguments[['params']] <- c(.self$fit_arguments[['params']],
                                             grid_of_hyperparameters)

      }





    },

    fit = function(data, ...) {

      x = sparse.model.matrix(formula(.self$formula),
                              data)

      if( .self$cross_validation){


      }else{


        xgb.mx <- xgb.DMatrix(data = x,
                              label = data[['deltaij']])

        setinfo(xgb.mx,
                        "base_margin",
                        log(data[['tij']]))

        .self$fit_arguments[['data']] <-xgb.mx

        out <- do.call(xgb.train,
                       .self$fit_arguments)

        }


      return(out)

    },
    predictor = function(model, newdata, ...) {


      # browser()
      x = sparse.model.matrix(formula(.self$formula),
                              newdata)

      # xgb.mx <- xgb.DMatrix(data = x,
      #                       label = newdata[['deltaij']])

      xgtest = xgb.DMatrix(x)
      setinfo(xgtest, "base_margin", log(newdata[['tij']]))

      out <- predict(model,
                     newdata = xgtest)

      return(out)


    },
    private_predictor = function(model, newdata, ...) {


      x = sparse.model.matrix(formula(.self$formula),
                              newdata)

      # xgb.mx <- xgb.DMatrix(data = x,
      #                       label = newdata[['deltaij']])

      xgtest = xgb.DMatrix(x)
      setinfo(xgtest, "base_margin", rep(1,nrow(newdata)))

      out <- predict(model,
                     newdata = xgtest)

      return(out)


    }




    ))

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

    fit = function(data, ...) {


      out <- glm.fit(x=sparse.model.matrix(formula(.self$formula),
                                           data),#[ss,]
                     y=data[['deltaij']], #[ss]
                     offset = log(data[['tij']]), #[ss]
                     family = poisson())

      return(out)

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


    predictor = function(model, newdata,...) {

      newdata<-newdata[complete.cases(newdata),]

      newx=sparse.model.matrix(formula(.self$formula),
                               newdata)

      newoffset = log(newdata[['tij']])

      out <-exp(newx %*% model$coefficients + newoffset)


      return(as.numeric(out))



    },

    private_predictor = function(model, newdata, ...) {



      newdata<-newdata[complete.cases(newdata),]


      newx=sparse.model.matrix(formula(.self$formula),
                               newdata)#[ss,]

      newoffset = log(newdata[['tij']])#[ss])


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




    update_cross_validation_argument= function(nfold){

      .self$fit_arguments[['nfolds']] <- nfold

    },

    fit = function(data, ...) {


      data<-data[complete.cases(data),]

      x = sparse.model.matrix(formula(.self$formula),
                              data)

      .self$fit_arguments[['y']] <- as.numeric(data[['deltaij']])

      .self$fit_arguments[['offset']] <-  log(data[['tij']])


      if(!.self$penalise_nodes){

        .self$fit_arguments[['penalty.factor']] <- 1- (grepl("node",colnames(x))& !grepl("node:", colnames(x))& !grepl(":node", colnames(x)))

      }

      .self$fit_arguments[['x']] <- x

      out <- do.call(.self$learner,
                     .self$fit_arguments)

      return(out)

    },


    predictor = function(model, newdata, ...) {


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
                     newoffset = log(1),#log(newdata[['tij']]),#[ss]),
                     type = "response")

      return(out)


    }



  )
)

#' \code{hal} learner class
#'
#' @export Learner_hal
#' @exportClass Learner_hal
Learner_hal <- setRefClass(
  "Learner_hal",
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
      .self$formula <- create_formula_hal(covariates = .self$covariates,
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



    update_cross_validation_argument= function(nfold){

      .self$fit_arguments[['nfolds']] <- nfold

    },

    fit = function(data, ...) {

      data <- data[, .(tij = sum(tij), deltaij = sum(deltaij)), by = c(c(covariates,treatment), "node", "k")]

      data<-data[complete.cases(data),]

      x = sparse.model.matrix(formula(.self$formula),
                              data)

      .self$fit_arguments[['y']] <- as.numeric(data[['deltaij']])

      .self$fit_arguments[['offset']] <-  log(data[['tij']])


      if(!.self$penalise_nodes){

        .self$fit_arguments[['penalty.factor']] <- 1- (grepl("node",colnames(x))& !grepl("node:", colnames(x))& !grepl(":node", colnames(x)))

      }

      .self$fit_arguments[['x']] <- x

      out <- do.call(.self$learner,
                     .self$fit_arguments)

      return(out)

    },


    predictor = function(model, newdata, ...) {



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
                     newoffset = log(1),#log(newdata[['tij']]),#[ss]),
                     type = "response")

      return(out)


    }



  )
)


#' \code{gam} learner class using mgcv bam
#'
#' @export Learner_gam
#' @exportClass Learner_gam
Learner_gam <- setRefClass(
  "Learner_gam",
  fields = list(
    covariates = "character",
    treatment = "character",
    cross_validation = "logical",
    recycle_information = "logical",
    intercept = "logical",
    formula = "character",
    learner = "function",
    add_nodes = "logical",
    penalise_nodes = "logical",
    fit_arguments = "list",
    covariates_attributes_matrix = "list"
  ),
  methods = list(

    initialize = function(covariates = NULL,
                          treatment = NA_character_,
                          cross_validation = FALSE,
                          intercept = TRUE,
                          add_nodes = TRUE,
                          penalise_nodes = FALSE,
                          recycle_information = FALSE,
                          ...) {

      .self$covariates <- covariates
      .self$treatment <- treatment
      .self$cross_validation <- cross_validation
      .self$intercept <- intercept
      .self$add_nodes <- add_nodes
      .self$penalise_nodes <- penalise_nodes
      .self$recycle_information <- recycle_information

      .self$formula <- create_formula_gam(
        covariates = .self$covariates,
        treatment = .self$treatment,
        intercept = .self$intercept,
        add_nodes = .self$add_nodes
      )

      .self$learner <- mgcv::bam

      .self$fit_arguments <- list(...)
      .self$fit_arguments[['family']] <- poisson()
    },

    datapp = function(train_data = NULL, validation_data = NULL) {
      if (!is.null(train_data)) {
        train.mf <- model.frame(as.formula(.self$formula),
                                rbind(train_data, validation_data),
                                drop.unused.levels = FALSE)
        x <- model.matrix(attr(train.mf, "terms"), data = train_data)
        y <- train_data[['deltaij']]
        offset <- log(train_data[['tij']])
        .self$covariates_attributes_matrix[['train.mf']] <- train.mf
        .self$recycle_information <- TRUE
      } else {
        x <- model.matrix(attr(.self$covariates_attributes_matrix[['train.mf']], "terms"),
                          data = validation_data)
        y <- validation_data[['deltaij']]
        offset <- log(validation_data[['tij']])
      }

      return(list(x = x, y = y, offset = offset))
    },

    fit = function(data, ...) {

      data <- data[complete.cases(data), ]
      .self$fit_arguments$formula <- as.formula(.self$formula)
      .self$fit_arguments$data <- data
      .self$fit_arguments$offset <- log(data[['tij']])
      fit <- do.call(.self$learner, .self$fit_arguments)
      return(fit)
    },


    predictor = function(model, newdata, ...) {

      pred <- predict(model,
                      newdata = newdata,
                      type = "response",
                      offset = log(newdata[['tij']]),
                      ...)
      return(pred)
    },

    private_predictor = function(model, newdata, ...) {


      pred <- predict(model,
                      newdata = newdata[node %in% model$xlevels$node,],
                      type = "response",
                      offset = log(1),#log(newdata[['tij']]),
                      ...)


      if(all((levels(newdata$node) %in% model$xlevels$node))){

        return(pred)

      }else{


        #
        newdata[node %in% model$xlevels$node,predictions_model:=pred]

        return(as.array(newdata$predictions_model))




        }
    }
  )
)
