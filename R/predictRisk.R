#' Absolute Risk prediction (Poisson Superlearner)
#'
#' This method computes the predicted absolute risk for some \code{cause}, at given \code{times} based on a \code{poisson_superlearner} object.
#'
#' @param object \code{poisson_superlearner} for absolute risk predictions.
#' @param newdata \code{data.frame}, new data to predict the absolute risk for.
#' @param times \code{numeric}, time(s) at which to predict the absolute risk.
#' @param cause \code{numeric}, competing risk to predict the absolute risk for.
#'
#' @return \code{matrix} with as many rows as the observations in \code{newdata} and as many columns as \code{times}.
#'
#'
#' @importFrom riskRegression predictRisk
#'
#' @references
#' Benichou, J., & Gail, M. H. (1990). Estimates of absolute cause-specific risk in cohort studies.
#' Biometrics, 46(3), 813–826. https://doi.org/10.2307/2532151 (PMID: 2242416)

#' @export
#' @export predictRisk.poisson_superlearner
predictRisk.poisson_superlearner <- function(object,newdata,times,cause=1, ...){

    rr_output <- matrix(NA_real_,
                        nrow=nrow(newdata),
                        ncol=length(times))


  for(ix in seq_along(times)){


    rr_output[,ix]<- predict(object,
            newdata = newdata,
            times = times[ix],
            cause = cause)$absolute_risk

  }

  # rr_output <- sapply(times, function(x){
  #     out <- predict(object,
  #             newdata = newdata,
  #             times = x,
  #             cause = cause)
  #
  #     p <- out$absolute_risk
  #
  #     return(p)
  #
  #   })

    # rr_output <- matrix(rr_output, nrow(newdata), ncol=length(times))


    return(rr_output)


}



#' Absolute Risk prediction (Learner)
#'
#' This method computes the predicted absolute risk for some \code{cause}, at given \code{times} based on a \code{base_learner} object.
#'
#' @param object \code{poisson_superlearner} for absolute risk predictions.
#' @param newdata \code{data.frame}, new data to predict the absolute risk for.
#' @param times \code{numeric}, time(s) at which to predict the absolute risk.
#' @param cause \code{numeric}, competing risk to predict the absolute risk for.
#'
#' @return \code{matrix} with as many rows as the observations in \code{newdata} and as many columns as \code{times}.
#'
#'
#' @importFrom riskRegression predictRisk
#'
#' @references
#' Benichou, J., & Gail, M. H. (1990). Estimates of absolute cause-specific risk in cohort studies.
#' Biometrics, 46(3), 813–826. https://doi.org/10.2307/2532151 (PMID: 2242416)

#' @export
#' @export predictRisk.poisson_superlearner
predictRisk.base_learner <- function(object,newdata,times,cause=1, ...){

  rr_output <- matrix(NA_real_,
                      nrow=nrow(newdata),
                      ncol=length(times))


  for(ix in seq_along(times)){


    rr_output[,ix]<- predict(object,
                             newdata = newdata,
                             times = times[ix],
                             cause = cause)$absolute_risk

  }

  # rr_output <- sapply(times, function(x){
  #     out <- predict(object,
  #             newdata = newdata,
  #             times = x,
  #             cause = cause)
  #
  #     p <- out$absolute_risk
  #
  #     return(p)
  #
  #   })

  # rr_output <- matrix(rr_output, nrow(newdata), ncol=length(times))


  return(rr_output)


}

