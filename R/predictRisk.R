#' Individual Data Pre-Processing
#'
#' Extract event probabilities from fitted poisson_superlearner.
#'
#' @return predictions
#'
#' @export
predictRisk.poisson_superlearner <- function(object,newdata,times,cause=1, ...){


    # browser()
    # rr_output <- matrix(nrow=nrow(newdata),
    #                     ncol=length(times))


    rr_output <- sapply(times, function(x){
      out <- predict(object,
              newdata = newdata,
              times = x,
              cause = cause)

      p <- 1 - out$survival_function

      return(p)

    })

    # out <- predict(object,
    #                newdata=newdata,
    #                times=times,
    #                cause=cause)

    # p <- 1 - out$survival_function
    #
    # return(p)

    rr_output <- matrix(rr_output, nrow(newdata), ncol=length(times))

    return(rr_output)


}
