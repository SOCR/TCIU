#' @title homologousX_features
#' @description Check if predictors are homologous
#' @details A general function that ensures the XReg predictors for 2 countries are homologous
#'
#' @param X_Country1 value inhereted from previous function "preprocess_ARIMA"
#' @param X_Country2 value inhereted from previous function "preprocess_ARIMA"
#'
#' @return this function returns a list of two data frames:
#' \item{X_Country1} {data frame with only X features from Country1. It shares the same features with Country2}
#' \item{X_Country2} {data frame with only X features from Country2. It shares the same features with Country1}
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#'
#' @export
#'
#' @examples
#' preprocess_Belgium <- preprocess_ARIMA( country = 'Belgium',
#' start=2000, end=2017, frequency=20,
#' feature="Gross domestic product at market prices")
#' preprocess_Bulgaria <- preprocess_ARIMA(country = 'Bulgaria',
#' start=2000, end=2017, frequency=20,
#' feature="Gross domestic product at market prices")
#'
#' homoFeat <- homologousX_features(preprocess_Belgium$X, preprocess_Bulgaria$X)
#' # See the two homologous data frames created by this function
#' summary(homoFeat)
#'
homologousX_features <- function (X_Country1, X_Country2){
  common_cols <- intersect(colnames(X_Country1), colnames(X_Country2))
  X_Country1 <- subset(X_Country1, select = common_cols)
  X_Country2 <- subset(X_Country2, select = common_cols)
  print(paste0("dim(X1)=(", dim(X_Country1)[1], ",", dim(X_Country1)[2], ");  ",
               " dim(X2)=(", dim(X_Country2)[1], ",", dim(X_Country2)[2], ")!"))
  return(list("X_Country1"=X_Country1, "X_Country2"=X_Country2))
}
