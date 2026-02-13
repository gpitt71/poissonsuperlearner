### test-glmnet-learners.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 12 2026 (06:30)
## Version:
## Last-Updated: feb 12 2026 (15:24)
##           By: Thomas Alexander Gerds
##     Update #: 18
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
library(testthat)
library(tmlensemble)
## library(glmnet)
library(survival)
library(riskRegression)
test_that("rigde works", {

    Xvars <- paste0("X", 1:10)
    d <- sampleData(n = 3000,formula = ~ f(X1, 2) + f(X2, 0) + f(X3, 0) + f(X6, 0) + f(X7, 0) + f(X8, 0) + f(X9, 0))

    fit_cox <- coxph(Surv(time, event == 1) ~ X1+X7,data = d,x = TRUE,y = TRUE)
    fit_pen_cox <- GLMnet(formula = Surv(time, event == 1) ~ X1+X7,data = d,lambda = 0)
    fit_rigde_cox <- GLMnet(formula = Surv(time, event == 1) ~ X1+X7,data = d,alpha = 0)
    fit_lasso_cox <- GLMnet(formula = Surv(time, event == 1) ~ X1+X7,data = d,alpha = 1)

    # Here you define the learner
    lridge <- Learner_glmnet(covariates = c("X1","X7"),
                              cross_validation=TRUE,
                              lambda_grid =seq(0.01,.9,by=0.002),
                              alpha=0,
                              intercept=TRUE,
                              penalise_nodes=FALSE)

    # Here you call the method that fits the learner to the data.

    # In my package you need explicit coding of the event.
    dupdated <- copy(d)
    dupdated[,event:=as.numeric(event==1)]

    lridge$fit(data = dupdated,
               event_time = "time",
               status = "event",
               number_of_nodes=20)

    # The model fit is saved as an attribute (model_fit) of the reference class Learner_glmnet.
    # It essentially works like a Python class.

    head(coef(lridge$model_fit[[1]]))


    lcox <- Learner_glmnet(covariates = c("X1","X7"),
                             cross_validation=FALSE,
                             lambda=0,
                             alpha=0,
                             intercept=TRUE,
                             penalise_nodes=FALSE)

    lcox$fit(data = dupdated,
               event_time = "time",
               status = "event",
             number_of_nodes = 20)

    head(coef(lcox$model_fit[[1]]))


    llasso <- Learner_glmnet(covariates = c("X1","X7"),
                             cross_validation=TRUE,
                             lambda_grid =seq(0.01,.9,by=0.002),
                             alpha=1,
                             intercept=TRUE,
                             penalise_nodes=FALSE)

    llasso$fit(data = dupdated,
             event_time = "time",
             status = "event",
             number_of_nodes = 20)



    lgam <- Learner_gam(covariates = c("X1","X7"))

    out_gam <- Superlearner(
      data = dupdated,
      event_time = "time",
      status = "event",
      id = "id",
      learners = list(lgam),
      nfold = 20,
      number_of_nodes = 20
    )

   cbind( lasso=head(coef(llasso$model_fit[[1]])),
    ridge=head(coef(lridge$model_fit[[1]])),
    cox=head(coef(lcox$model_fit[[1]])),
    gam=head(coef(out_gam$superlearner[[1]]$learners_fit)))


})


test_that("rigde works in Steno1", {
    covariates <- c("sex","age","diabetes_duration","value_SBP","value_LDL","value_HBA1C","value_Smoking","value_Motion","value_Albuminuria","eGFR")
    d <- simulateStenoT1(n = 300,beta_age_rate_cvd = .1,beta_age_rate_death = 0,keep = c("time","event",covariates))
    fit_cox <- coxph(Surv(time,event == 1)~age,data = d,x = TRUE,y = TRUE)
    lasso_fish <- fishNet(data = d,event_time = "time",status = "event",covariates = covariates,alpha = 1,lambda = seq(.001,.9,.0002),penalise_nodes = FALSE,number_of_nodes=20,nfold = 20)
    ridge_fish <- fishNet(data = d,event_time = "time",status = "event",covariates = covariates,alpha = 0,lambda = seq(.001,.9,.0002),penalise_nodes = FALSE,number_of_nodes=20,nfold = 20)
    np_ridge_fish <- fishNet(data = d,event_time = "time",status = "event",covariates = covariates,alpha = 0,lambda = seq(.001,.9,.0002),penalise_nodes = TRUE,number_of_nodes=20,nfold = 20)
    cbind(lasso_fish$superlearner[[1]]$learners_fit$beta,
          ridge_fish$superlearner[[1]]$learners_fit$beta,
          np_ridge_fish$superlearner[[1]]$learners_fit$beta)
})




test_that("glmnet convergence issues veteran", {

  veteran_data <- as.data.table(survival::veteran)
  veteran_data <- veteran_data[complete.cases(veteran_data), ]
  veteran_data[, id := 1:.N]
  veteran_data[,prior:=as.factor(prior)]

  veteran_data <- veteran_data[,.SD[1],by="time"] # disregard ties

  covariates <- c("age","celltype")

  fit_cox <- coxph(Surv(time,status )~age+celltype,data = veteran_data,x = TRUE,y = TRUE)

  fit_rigde_cox <- GLMnet(formula = Surv(time, status == 1) ~ age+celltype,data = veteran_data,alpha = 0)
  fit_lasso_cox <- GLMnet(formula = Surv(time, status == 1) ~ age+celltype,data = veteran_data,alpha = 1)

  lglmnet <- Learner_glmnet(covariates = covariates,
                                  cross_validation=FALSE,
                                  lambda=0,
                                  lambda_grid =seq(.001,.9,.0002),
                                  intercept=FALSE,
                                  penalise_nodes=FALSE)

  lglmnet$fit(veteran_data,
              id="id",
              status="status",
              stratified_k_fold=FALSE,
              event_time = "time")

  lglmnet_yt <- Learner_glmnet(covariates = covariates,
                            cross_validation=FALSE,
                            alpha=1,
                            lambda=0,
                            lambda_grid =seq(.001,.9,.0002),
                            intercept=FALSE,
                            thresh = 1e-15,
                            penalise_nodes=FALSE)


  lglmnet_yt$fit(veteran_data,
              id="id",
              status="status",
              stratified_k_fold=FALSE,
              event_time = "time")


  llasso <- Learner_glmnet(covariates = covariates,
                               cross_validation=T,
                               alpha=1,
                               # lambda=0,
                               lambda_grid =seq(.001,.9,.0002),
                               intercept=FALSE,
                               thresh = 1e-15,
                               penalise_nodes=TRUE)


  llasso$fit(veteran_data,
                 id="id",
                 status="status",
                 stratified_k_fold=FALSE,
                 event_time = "time")


  lgam <- Learner_gam(covariates = c("age","celltype"))

  out_gam <- Superlearner(veteran_data,
                                      id="id",
                                      status="status",
                                      stratified_k_fold=FALSE,
                                      event_time = "time",
                                      learners=list(lgam),
                                      meta_learner_algorithms = c("glm"),
                                      nfold = 20
  )

  cbind(coef(fit_cox),
        lglmnet$model_fit$beta[1:5,],
        lglmnet_yt$model_fit$beta[1:5,],
        llasso$model_fit$beta[1:5,],
        out_gam$superlearner[[1]]$learners_fit$coefficients[2:5]
  )

})

######################################################################
### test-glmnet-learners.R ends here
