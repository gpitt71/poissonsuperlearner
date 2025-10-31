#' Absolute Risk prediction
#'
#' This method computes the predicted absolute risk for some \code{cause}, at given \code{times} based on a \code{poisson_superlearner} object.
#'
#' @param object \code{poisson_superlearner} for absolute risk predictions.
#' @param newdata \code{data.frame}, new data to predict the absolute risk for.
#' @param times \code{numeric}, time(s) at which to predict the absolute risk.
#' @param cause \code{numeric}, competing risk to predict the absolute risk for.
#' @param absolute_risk_integration \code{character}. Using the argument \code{"exact"}, the function computes the predicted absolute risk as in Benichou and Gail (1990). The argument \code{"approx"} computes Equation 1.1. of Benichou and Gail (1990) as a discrete sum.
#'
#' @return \code{matrix} with as many rows as the observations in \code{newdata} and as many columns as \code{times}.
#'
#' @references
#' Benichou, J., & Gail, M. H. (1990). Estimates of absolute cause-specific risk in cohort studies.
#' Biometrics, 46(3), 813–826. https://doi.org/10.2307/2532151 (PMID: 2242416)

#' @export
predictRisk.poisson_superlearner <- function(object,newdata,times,cause=1,absolute_risk_integration="exact", ...){

    rr_output <- matrix(NA_real_,
                        nrow=nrow(newdata),
                        ncol=length(times))


  for(ix in seq_along(times)){


    rr_output[,ix]<- predict(object,
            newdata = newdata,
            times = times[ix],
            absolute_risk_integration=absolute_risk_integration,
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
