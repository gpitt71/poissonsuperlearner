### print.fishNet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 12 2026 (08:30) 
## Version: 
## Last-Updated: feb 12 2026 (08:46) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' print fishNet objects
##'
##' show beta coefficients
##' @title print fishNet objects
##' @param x object obtained with \link[tmlensemble]{fishNet}.
##' @param ... passed to print
##' @return object
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
print.fishNet <- function(x,...){
    if (length(x$superlearner[[1]]$learners_fit$glmnet.fit)>0){
        print(x$superlearner[[1]]$learners_fit$glmnet.fit$beta)
    }else{
        print(x$superlearner[[1]]$learners_fit$beta)
    }
    invisible(x)
}


######################################################################
### print.fishNet.R ends here
