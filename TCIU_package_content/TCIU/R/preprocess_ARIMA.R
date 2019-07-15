#' @title preprocess_ARIMA
#' @description  Preprocess dataset for ARIMA model
#' @details This function preprocesses dataset of given countries to ensure full rank
#' Then, extend the Fit-ARIMA method to ensure testing-training
#' modeling/assessment for 2 countries works
#'
#' @param country Give a country name that is shown on the original dataset
#' @param start Select the start year to create the ARIMA model
#' @param end Select the end year to create the ARIMA model
#' @param frequency The number of observations per unit of time. The same in function "ts"
#' @param feature Choose one feature to create the ARIMA model on
#' 
#' @return a list contains two elements:
#' \item{X}{A data frame contains the full rank part of all features}
#' \item{Y}{A vector contains response variable}
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
#' 
#' @examples 
#' preprocess_Belgium <- preprocess_ARIMA( country = 'Belgium', 
#' start=2000, end=2017, frequency=20, 
#' feature="Gross domestic product at market prices")
#' 
#' # Just show the response variable part
#' head(preprocess_Belgium$Y)
preprocess_ARIMA <-
  function(country='Belgium',start=2000, end=2017, frequency=20,
           feature="Unemployment , Females, From 15-64 years, Total"){
    # delete features that are missing for all time points
    DataSuperSample = rmv_miss_ftr(countryName = country)

    print(paste0("Processing feature: ...", feature, "... "))

    if (feature == "Unemployment , Females, From 15-64 years, Total") {
      Y = select(DataSuperSample,
                 "Unemployment , Females, From 15-64 years, Total")
      X = select(DataSuperSample,
                 -starts_with("Unemployment"))
    } else if (feature == "Gross domestic product at market prices") {
      Y = select(DataSuperSample,
                 "Gross domestic product at market prices"); dim(Y)
      X = select(DataSuperSample,
                 -matches("debt|Debt") )
      X <- X [, -c(50:80)]; dim(X)  # 360 167
    } else {
      print(paste0("This feature: ...",
                   feature,
                   "... is not imlemented yet! Exiting preprocess_ARIMA() method ..."))
      return(NULL)
    }
    # reduce the number of observations (matrix rows) to specified timerange
    len_1 <- (end + 1 - start) * frequency; print(paste0("dim(X)[1]=", len_1))
    X <- X[1:len_1 , qr(X[1:len_1 , ])$pivot[seq_len(qr(X[1:len_1 , ])$rank)]]; dim(X)
    # ensure full-rank design matrix, X
    Y <- as.data.frame(Y[1:len_1 , ])
    print(paste0("dim(X)=(", dim(X)[1], ",", dim(X)[2], ");  ",         # 300 136
                 " dim(Y)=(", dim(Y)[1], ",", dim(Y)[2], ") ..."))    # 300 1
    return( list("X"=X, "Y"=Y) )
  }
