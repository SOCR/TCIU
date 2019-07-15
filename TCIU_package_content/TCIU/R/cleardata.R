#' @title cleardata
#' @description This function imputes missing values for a dataset
#'
#' @param mat A matrix to be cleaned.
#'
#' @details This function use the mean value of each column to impute each missing value (i.e. NA value)
#' in that column. Thus provides a quick and easy method to resolve value missing problem in a dataset.
#'
#' @return A matrix whose NA's are all filled with column mean of the rest items.
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#'
#' @examples
#' #Need package "VIM" to make this example work
#' require(VIM)
#'
#' # See the missing value in dataset "sleep"
#' head(VIM:sleep)
#'
#' # Deal with missing values
#' cleardata(VIM:sleep)
#'
#' @export
#'
cleardata <- function(mat) {
  for (i in 1:ncol(mat)) {
    mat[is.na(mat[,i]),i]<-mean(mat[,i],na.rm = T)
  }
  #remove data with only NA/NAN value
  mat = mat[,which(is.nan(mat[1,])==FALSE)]

  ## Below is a version that utilizes R's vectorization property
  # mat =
  #   apply(mat, 2, function(column)
  #     sapply(column, function(x)
  #       ifelse(is.na(x),
  #              mean(column,na.rm = T) +
  #                rnorm( sum(is.na(column) ), sd = sd(column, na.rm = T) ),
  #              x) ) )
  return(mat)
}
