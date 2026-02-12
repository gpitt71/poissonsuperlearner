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
    Xvars <- paste("X",1:10)
    d <- sampleData(n = 300,formula = ~f(X1, 2) + f(X2, 0) + f(X3, 0) + f(X6, 0) + f(X7, 0) + f(X8, 0) + f(X9, 0))
    fit_cox <- coxph(Surv(time,event == 1)~X1,data = d,x = TRUE,y = TRUE)
    lasso_fish <- fishNet(data = d,event_time = "time",status = "event",covariates = Xvars,alpha = 1,lambda = seq(.001,.9,.0002),penalise_nodes = FALSE,number_of_nodes=20,nfold = 20)
    ridge_fish <- fishNet(data = d,event_time = "time",status = "event",covariates = Xvars,alpha = 0,lambda = seq(.001,.9,.0002),penalise_nodes = FALSE,number_of_nodes=20,nfold = 20)
    cbind(lasso_fish$superlearner[[1]]$learners_fit$beta,
          ridge_fish$superlearner[[1]]$learners_fit$beta)
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


######################################################################
### test-glmnet-learners.R ends here
