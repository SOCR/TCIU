#' @title noiseGen
#' @description This function is useful for generating noises for a given dataset
#' @details Noise-generating function to generate some synthetic predictors (X_new)
#' and test the ARIMA(2,0,1)(1,0,1)[12] model
#'
#'
#'
#' @param data A data in the form of vector or of a matrix
#' @return A updated version of data where noise is also included
#' @examples
#' noiseGen( matrix(rep(2,9),3,3) )
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
noiseGen <- function(data) {
  if (is.vector(data))
    return( data + rnorm(length(data), mean(data), sd(data)))
  else
    return(data + matrix(rnorm(dim(data)[1] * dim(data)[2], 0, 1), dim(data)[1]))
}


