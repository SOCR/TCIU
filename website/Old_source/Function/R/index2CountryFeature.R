#' @title index2CountryFeature
#' @description Convert index to country name and feature code
#' @details This function maps to convert between 1D indices
#' \strong{Notice:} This function mainly serves as a functional part for other functions. So we seldomly
#' apply this function by itself.
#'
#' @param indx Indicate the feature you wish to see after transfering to 1D index
#' @return return the original country and feature code of the feature you wish to show
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#'
#' @export
#'
index2CountryFeature <- function(indx=1) {
  if (indx<1 | indx>length(arimaModels_ARMA_coefs)) {
    cat("Index out of bounds: indx=", indx, "; must be between 1 and ",
        length(arimaModels_ARMA_coefs), " ... Exiting ...")
    return (NULL)
  } else {
    feature = (indx-1) %% (dim(list_of_dfs_CommonFeatures[[1]])[2])
    country = floor((indx - feature)/(dim(list_of_dfs_CommonFeatures[[1]])[2]))
    return(list("feature"=(feature+1), "country"=(country+1)))  }
}
