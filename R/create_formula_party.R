#' Build a glmnet Poisson formula
#'
#' Developer note: internal formula constructor for the `Learner_glmnet`
#' design matrix path.
#'
#' Assumes covariate names are syntactically valid formula terms. The node term
#' is always added and the exposure offset is included. Invalid names or missing
#' required columns fail downstream during model-matrix construction.
#'
#' @param covariates Character vector of covariate names or `NA_character_`.
#' @returns Character scalar formula for `deltaij`.
#' @noRd
create_formula_glmnet <- function(covariates=NA_character_){


  xs<- NULL

  if (!any(is.na( covariates))) {
    xs <- paste(covariates, collapse = "+")
  }

  xs <- paste(xs, "+ node")


  out <- paste("deltaij ~", xs, "+offset(log(tij))", sep =
                 "")

  return(out)




}

#' Build a GAM formula
#'
#' Developer note: internal formula constructor for `Learner_gam`.
#'
#' Assumes covariate names are valid formula terms. The node term is always
#' added, and no exposure offset is included here because the backend handles it
#' separately. The `competing_risks` argument is retained for compatibility.
#'
#' @param covariates Character vector of covariate names or `NA_character_`.
#' @param competing_risks Historical flag, currently unused.
#' @returns Character scalar formula for `deltaij`.
#' @noRd
create_formula_gam <- function(covariates=NA_character_,
                           competing_risks=FALSE){



  if (!any(is.na( covariates))) {
    xs <- paste(covariates, collapse = "+")
  }

  xs <- paste(xs, "+ node")


  out <- paste("deltaij ~", xs, sep =
                 "")

  return(out)




}

#' Build a HAL interaction formula
#'
#' Developer note: internal formula constructor for HAL-style basis generation.
#'
#' Assumes covariate names are valid formula terms. Covariates and `node` are
#' joined with `*` so interaction terms can be expanded by downstream tooling.
#' The `competing_risks` argument is retained for compatibility. Invalid names
#' fail later during formula parsing or basis construction.
#'
#' @param covariates Character vector of covariate names or `NA_character_`.
#' @param competing_risks Historical flag, currently unused.
#' @param intercept Logical; whether to keep an intercept.
#' @returns Character scalar formula with exposure offset.
#' @noRd
create_formula_hal <- function(covariates=NA_character_,
                           competing_risks=FALSE,
                           intercept=FALSE){



  if (!any(is.na( covariates))) {
    xs <- paste(covariates, collapse = "*")
  }

  # if (!is.na( treatment)) {
  #   xs <- paste(xs, "*", treatment)
  # }

  # if (competing_risks) {
  #   xs <- paste(xs, "+ k")
  # }

  xs <- paste(c(xs, "node"), collapse = "*")


  if(!intercept){
    xs <- paste(xs, "-1")
  }

  out <- paste("deltaij ~", xs, "+offset(log(tij))", sep =
                 "")

  return(out)




}
