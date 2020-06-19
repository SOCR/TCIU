#' @title coef_plot
#' @description Plotting the coefficients
#' @details This function plots the coefficients
#' \strong{Notice:} This function mainly serves as a functional part for other functions. So we seldomly
#' apply this function by itself.
#'
#' @param betahat A vector that stores estimated coefficients.
#' @param varn A vector that stores names for the coefficients.
#' @param plotname A string that to be used as the plot name.
#' @return A plot of coefficients.
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#'
#' @export
#'
coef_plot <- function(betahat,varn,plotname) {
  betahat <- betahat[-1]
  P <-
    coefplot(betahat[which(betahat!=0)],
             sd = rep(0, length(betahat[which(betahat!=0)])),
             pch=0, cex.pts = 3, col="red", main = plotname,
             varnames = varn[which(betahat!=0)] )
  return(P)
}
