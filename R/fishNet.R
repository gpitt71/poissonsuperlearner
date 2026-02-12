### fishNet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 12 2026 (08:18) 
## Version: 
## Last-Updated: feb 12 2026 (08:33) 
##           By: Thomas Alexander Gerds
##     Update #: 5
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Fit poisson regression model  
##'
##' via Superlearner
##' @title fit poisson regression model
##' @param data data including the outcome variables and the covariates
##' @param event_time character. event time variable 
##' @param status character. event type. 
##' @param covariates character vector of covariate names
##' @param alpha alpha mixing parameter (1=LASSO, 0=ridge, elastic net else)
##' @param lambda_grid vector of values for lambda search
##' @param penalise_nodes if TRUE penalize also baseline hazard parameters 
##' @param number_of_nodes number of baseline hazard parameters
##' @param nfold number of folds for super learner
##' @return fitted object
##' @examples
##' covariates <- c("sex","age","diabetes_duration","value_SBP","value_LDL",
##'     "value_HBA1C","value_Smoking","value_Motion","value_Albuminuria","eGFR")
##' d <- simulateStenoT1(n = 300,beta_age_rate_cvd = .1,beta_age_rate_death = 0,
##'                      keep = c("time","event",covariates))
##' fit <- fishNet(data = d,event_time = "time",status = "event",
##'  covariates = covariates,alpha = 0,lambda = seq(.001,.9,.0002),
##'  penalise_nodes = FALSE,number_of_nodes=20,nfold = 20)
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
fishNet <- function(data,
                    event_time,
                    status,
                    covariates,
                    alpha,
                    lambda_grid,
                    penalise_nodes = FALSE,
                    number_of_nodes=20,
                    nfold = 20){
    fit <- Learner_glmnet(covariates,
                          cross_validation=TRUE,
                          alpha=1,
                          lambda_grid =lambda_grid,
                          intercept=FALSE)
    sl_fit <- suppressWarnings(Superlearner(data,
                                            stratified_k_fold=FALSE,
                                            event_time = event_time,
                                            status=status,
                                            learners=list(fit),
                                            meta_learner_algorithms = c("glm"),
                                            number_of_nodes = number_of_nodes,
                                            nfold = nfold))
    class(sl_fit) <- "fishNet"
    sl_fit
}
######################################################################
### fishNet.R ends here
