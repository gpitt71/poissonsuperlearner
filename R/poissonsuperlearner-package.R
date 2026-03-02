#' @keywords internal
"_PACKAGE"

#' @importFrom riskRegression predictRisk
#' @importFrom methods new
#' @importFrom stats coef deviance predict quantile relevel time
#' @importFrom utils head tail
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom mgcv gam
#' @importFrom Matrix sparse.model.matrix
#' @import data.table
#' @import sampling
NULL
## usethis namespace: start
#' @useDynLib poissonsuperlearner, .registration = TRUE
## usethis namespace: end
NULL
