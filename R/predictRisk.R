#' Individual Data Pre-Processing
#'
#' Extract event probabilities from fitted poisson_superlearner.
#'
#' @return predictions
#'
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
