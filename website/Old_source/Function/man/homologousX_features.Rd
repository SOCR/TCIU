% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homologousX_features.R
\name{homologousX_features}
\alias{homologousX_features}
\title{homologousX_features}
\usage{
homologousX_features(X_Country1, X_Country2)
}
\arguments{
\item{X_Country1}{value inhereted from previous function "preprocess_ARIMA"}

\item{X_Country2}{value inhereted from previous function "preprocess_ARIMA"}
}
\value{
this function returns a list of two data frames:
\item{X_Country1} {data frame with only X features from Country1. It shares the same features with Country2}
\item{X_Country2} {data frame with only X features from Country2. It shares the same features with Country1}
}
\description{
Check if predictors are homologous
}
\details{
A general function that ensures the XReg predictors for 2 countries are homologous
}
\examples{
preprocess_Belgium <- preprocess_ARIMA( country = 'Belgium',
start=2000, end=2017, frequency=20,
feature="Gross domestic product at market prices")
preprocess_Bulgaria <- preprocess_ARIMA(country = 'Bulgaria',
start=2000, end=2017, frequency=20,
feature="Gross domestic product at market prices")

homoFeat <- homologousX_features(preprocess_Belgium$X, preprocess_Bulgaria$X)
# See the two homologous data frames created by this function
summary(homoFeat)

}
\author{
SOCR team <http://socr.umich.edu/people/>
}
