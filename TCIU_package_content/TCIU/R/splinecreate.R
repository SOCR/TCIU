#' @title splinecreate
#' @description Expand the dataset using spline regression
#' 
#' @details  This function uses spline regression model to expand the dataset
#' Here we will use the spline model to create independent smooth regression model
#' for each features. Then estimate new prediction points inside the original time
#' series points. Thus expands the original quarterly data into a weekly data.
#' During the process, a Gaussian noise will be added onto prediction values to
#' make the prediction data more reasonable.
#' To deal with the problem of rank deficiency, we'll use full-rank matrix in the ARIMAX model.
#' This means some more features will be eliminated from the original model.
#' Moreover, a parameter k would be used to narrow down the variance of the
#' Gaussian noise we apply artifically, in order to get a better-fitted model.
#'
#' @param mat A matrix.
#' @param k A number that helps control the variance of the Gaussian noise.
#' @return A matrix that contains the regression spline results for each feature.
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
#' 
#' @examples 
#' # Plot with original dataset (also take the first column as an example)
#' plot(Belgium_Clear[,1],ylab="Value")
#' 
#' # Plot with same column but using gaussian noise to expand the dataset
#' # Take k=1 then we are using one standard error
#' head(splinecreate(Belgium_Clear,1))
#' plot(splinecreate(Belgium_Clear,1)[,1],ylab="Value")
splinecreate <- function(mat,k) {
  res<-NULL
  for (i in 1:ncol(mat)) {
    sp <- smooth.spline(seq(1:72),mat[,i])
    spresult <- predict(sp,seq(from=0,by=1/13,length.out = 936))
    spfeat <- spresult$y + rnorm(936,sd=(sd(spresult$y))/k)
    res <- cbind(res,spfeat)
  }
  colnames(res) <- colnames(mat)
  return(res)
}
