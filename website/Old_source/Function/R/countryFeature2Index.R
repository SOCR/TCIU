#' @title countryFeature2Index
#' @description Convert country feature to index
#' @details  This function transfer country and feature code to 1D index
#' An opposite procedure of function "index2CountryFeature"
#' \strong{Notice:} This function mainly serves as a functional part for other functions. So we seldomly
#' apply this function by itself.
#'
#' @param countryIndx indicate the country you wish to transfer
#' @param featureIndx indicate the feature you wish to transfer
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#'
#' @export
#'
countryFeature2Index <- function(countryIndx=1, featureIndx=1) {
  if (countryIndx<1 | countryIndx>(dim(list_of_dfs_CommonFeatures[[1]])[1]) |
      featureIndx<1 | featureIndx>(dim(list_of_dfs_CommonFeatures[[1]])[2])) {
    cat("Indices out of bounds: countryIndx=", countryIndx, "; featureIndx=", featureIndx, "; Exiting ...")
    return (NULL)
  } else { return (featureIndx + (countryIndx-1)*(dim(list_of_dfs_CommonFeatures[[1]])[2])) }
}
