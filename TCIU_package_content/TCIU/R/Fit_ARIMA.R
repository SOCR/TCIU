#' @title Fit_ARIMA
#' @description  Fit the ARIMA model
#' @details This function cleans the data and fits the ARIMA model. It is based on several previous functions
#' and must work under the European Union Economy dataset. It is built primarily aiming to create ARIMA model for 
#' the EU Economy dataset.
#'
#' @param country Give a country name that is shown on the original dataset
#' @param start Select the start year to create the ARIMA model
#' @param end Select the end year to create the ARIMA model
#' @param frequency The number of observations per unit of time. The same in function "ts"
#' @param feature Choose one feature to create the ARIMA model on
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
#' 
#' @examples 
#' # Example are based on country 
#' BulgariaARIMA = 
#' Fit_ARIMA(country ='Bulgaria', 
#'          start=2000, end=2017, frequency=20,
#'          feature="Gross domestic product at market prices")
#' summary(BulgariaARIMA)
Fit_ARIMA <- function(country = 'Belgium',
                      start = 2000, end = 2017, frequency = 20,
                      feature="Unemployment , Females, From 15-64 years, Total")
{
  # delete features that are missing for all time points
  DataSuperSample = rmv_miss_ftr(countryName = country)
  if (feature=="Unemployment , Females, From 15-64 years, Total") {
    Y = select(DataSuperSample, "Unemployment , Females, From 15-64 years, Total")
    X = select(DataSuperSample, -starts_with("Unemployment"))
  } else if (feature=="Gross domestic product at market prices") {
    Y = select(DataSuperSample, "Gross domestic product at market prices"); dim(Y)
    X = select(DataSuperSample, -matches("debt|Debt")); dim(X)  # 360 167
    print(paste0("dim(X)=(", dim(X)[1], ",", dim(X)[2], ");  ",
                 " dim(Y)=(", dim(Y)[1], ",", dim(Y)[2], ") ..."))
    # ensure full-rank design matrix, X
    X <- X[, qr(X)$pivot[seq_len(qr(X)$rank)]]; dim(X)
  }
  else {
    print(paste0("This feature ", feature, " is not imlemented yet! Exiting Fit_ARIMA() method ..."))
    return(NULL)
  }
  ts_Y <- ts(Y, start=c(start, 1), end=c(end, frequency),
             frequency = frequency); length(ts_Y)
  set.seed(1234)
  fitArimaX = auto.arima(ts_Y, xreg=as.matrix(X) )
  return(fitArimaX)
}
