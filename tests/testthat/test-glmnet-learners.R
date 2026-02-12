### test-glmnet-learners.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 12 2026 (06:30)
## Version:
## Last-Updated: feb 12 2026 (09:05)
##           By: Thomas Alexander Gerds
##     Update #: 12
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
  d <- sampleData(
    n = 300,
    formula = ~ f(X1, 2) + f(X2, 0) + f(X3, 0) + f(X6, 0) + f(X7, 0) + f(X8, 0) + f(X9, 0)
  )
  fit_cox <- coxph(Surv(time, event == 1) ~ X1,
                   data = d,
                   x = TRUE,
                   y = TRUE)
  lasso_fish <- fishNet(
    data = d,
    event_time = "time",
    status = "event",
    covariates = Xvars,
    alpha = 1,
    lambda = seq(.001, .9, .0002),
    penalise_nodes = FALSE,
    number_of_nodes = 20,
    nfold = 20
  )
  ridge_fish <- fishNet(
    data = d,
    event_time = "time",
    status = "event",
    covariates = Xvars,
    alpha = 0,
    lambda = seq(.001, .9, .0002),
    penalise_nodes = FALSE,
    number_of_nodes = 20,
    nfold = 20
  )
  cbind(
    lasso_fish$superlearner[[1]]$learners_fit$beta,
    ridge_fish$superlearner[[1]]$learners_fit$beta
  )
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

  covariates <- c("age")

  fit_cox <- coxph(Surv(time,status )~age,data = veteran_data,x = TRUE,y = TRUE)

  lglmnet <- Learner_glmnet(covariates = c("age"),
                                  cross_validation=FALSE,
                                  alpha=1,
                                  lambda=0,
                                  lambda_grid =seq(.001,.9,.0002),
                                  intercept=FALSE,
                                  penalise_nodes=FALSE)

  out_glmnet_notresh <- Superlearner(veteran_data,
                                   id="id",
                                   status="status",
                                   stratified_k_fold=FALSE,
                                   event_time = "time",
                                   learners=list(lglmnet),
                                   meta_learner_algorithms = c("glm"),
                                   nfold = 20
  )

  lglmnet_yt <- Learner_glmnet(covariates = c("age"),
                            cross_validation=FALSE,
                            alpha=1,
                            lambda=0,
                            lambda_grid =seq(.001,.9,.0002),
                            intercept=FALSE,
                            thresh = 1e-15,
                            penalise_nodes=FALSE)

  out_glmnet_yestresh <- Superlearner(veteran_data,
                                     id="id",
                                     status="status",
                                     stratified_k_fold=FALSE,
                                     event_time = "time",
                                     learners=list(lglmnet_yt),
                                     meta_learner_algorithms = c("glm"),
                                     nfold = 20
  )

  lgam <- Learner_gam(covariates = c("age"))

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
  out_glmnet_notresh$superlearner[[1]]$learners_fit$beta[1,],
  out_glmnet_yestresh$superlearner[[1]]$learners_fit$beta[1,],
  out_gam$superlearner[[1]]$learners_fit$coefficients[2]
  )

})

######################################################################
### test-glmnet-learners.R ends here
