#' @title forecast the fMRI data based on the time series
#' @description a function to forecast the fMRI data based on the time series
#'
#' @param fmridata a 4D array contains information for the fMRI spacetime image. The data should only contain the magnitude for the fMRI image.
#' @param voxel_location a 3d array indicating the voxel location of the brain
#' @param cut breaking point of the time-series data. The default is 10. 
#'
#' @details The function \code{fmri_ts_forecast} is used to forecast with time series. It will fit the best ARIMA model to univariate time series from the input fMRI data.
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return a figure forecasting the fMRI voxel with time series
#' @export
#' @import zoo plotly
#' @importFrom stats time
#' @importFrom forecast auto.arima forecast
#' 
#' @examples
#' require('AnalyzeFMRI')
#' fmri_generate = fmri_simulate_func(dim_data = c(64, 64, 40), mask = mask)
#' \donttest{
#' smoothmod <- AnalyzeFMRI::GaussSmoothArray(fmri_generate$fmri_data, sigma = diag(3,3))
#' fmri_ts_forecast(smoothmod,c(41,44,33))
#' }

fmri_ts_forecast = function(fmridata,
                            voxel_location,
                            cut = 10) {
    data = detrend(fmridata[voxel_location[1], voxel_location[2], voxel_location[3], ], bp = seq(21, 160, by = cut))
    trainData <- data[3:132]
    testData <- ts(data[133:160], start = 131, end = 158)
    arimaMod <- auto.arima(trainData, stepwise = FALSE, approximation = FALSE)
    TSplot_gen(130, arimaMod, periods = 28, ts_list = list(testData))
}
