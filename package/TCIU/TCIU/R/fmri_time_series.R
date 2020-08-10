#' @title visualization of the fMRI data (real, imaginary, magnitude, and phase parts) in time series
#' @description a visualization method, use \code{plotly} to draw the fMRI data in time series
#'
#' @param fmridata a 4d array which contains the spatial and temporal record of fMRI result or a single complex valued vector
#' @param voxel_location a 3d array indicating the spatial location of the brain. If is.4d is false, set the voxel_location as NULL.
#' @param is.4d The default is TRUE. If change to false, input a vector instead of a 4d array.
#' @param ref The default is NULL. User can input an outside extra reference \code{plotly} object to include in the final result.
#'
#' @details
#' The function \code{fmri_time_series} is used to create four interactive time series graphs for the real, imaginary, magnitude, and phase parts for the fMRI spacetime data.
#' User can choose to provide the 4d array of the fMRI spacetime image and the voxel_location or a single complex valued vector, then four interactive time series graphs will be shown. 
#' Besides, the reference \code{plotly} object can be added to the final result.
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return an interactive time series graph object created by \code{plotly}
#' @export
#'
#' @import plotly
#' @importFrom pracma detrend
#' @examples
#' # load sample time-series data of one voxel in the brain provided by the package
#' sample_voxel = sample[[5]]
#' reference_plot = sample[[4]]
#' fmri_time_series(sample_voxel, voxel_location = NULL, is.4d = FALSE, ref = reference_plot)
#' 

fmri_time_series = function(fmridata,
                            voxel_location,
                            is.4d = TRUE,
                            ref = NULL) {
    if (is.4d) {
        X = voxel_location[1]
        Y = voxel_location[2]
        Z = voxel_location[3]
        realnum <- Re(fmridata[X, Y, Z, ])
        imgnum <- Im(fmridata[X, Y, Z, ])
        phasenum <- Arg(fmridata[X, Y, Z, ])
        modnum <- Mod(fmridata[X, Y, Z, ])
    } else {
        realnum <- Re(fmridata)
        imgnum <- Im(fmridata)
        phasenum <- Arg(fmridata)
        modnum <- Mod(fmridata)
    }
    realnum1 <- detrend(realnum, bp = seq(21, 160, by = 20))
    tsrealnum1 <- ts(realnum1)
    ksmthrealnum <- ksmooth(c(1:160), tsrealnum1, kernel = "normal", bandwidth = 5)
    ksthrealnum <- data.frame(tsrealnum1, ksmthrealnum$y)
    TScore_realnum <- GTSplot(ksthrealnum, Unit = "time point", ts_name = c("real_original", "real_ksmooth"), COLO = c("FFCC33", 
        "00CCFF"))
    imgnum1 <- detrend(imgnum, bp = seq(21, 160, by = 20))
    tsimgnum1 <- ts(imgnum1)
    ksmthimgnum <- ksmooth(c(1:160), tsimgnum1, kernel = "normal", bandwidth = 5)
    ksthimgnum <- data.frame(tsimgnum1, ksmthimgnum$y)
    TScore_imgnum <- GTSplot(ksthimgnum, Unit = "time point", ts_name = c("img_original", "img_ksmooth"), COLO = c("FF9966", "0099FF"))
    phasenum1 <- detrend(phasenum, bp = seq(21, 160, by = 20))
    tsphasenum1 <- ts(phasenum1)
    ksmthphasenum <- ksmooth(c(1:160), tsphasenum1, kernel = "normal", bandwidth = 5)
    ksthphasenum <- data.frame(tsphasenum1, ksmthphasenum$y)
    TScore_phasenum <- GTSplot(ksthphasenum, Unit = "time point", ts_name = c("phase_original", "phase_ksmooth"), COLO = c("FF6633", 
        "0066FF"))
    modnum1 <- detrend(modnum, bp = seq(21, 160, by = 20))
    tsmodnum1 <- ts(modnum1)
    ksmthmodnum <- ksmooth(c(1:160), tsmodnum1, kernel = "normal", bandwidth = 5)
    ksthmodnum <- data.frame(tsmodnum1, ksmthmodnum$y)
    TScore_modnum <- GTSplot(ksthmodnum, Unit = "time point", ts_name = c("mod_original", "mod_ksmooth"), COLO = c("CC3300", "0033FF"))
    if (is.null(ref)) {
        result <- plotly::subplot(TScore_realnum, TScore_imgnum, TScore_phasenum, TScore_modnum, nrows = 4)
    } else {
        result <- plotly::subplot(TScore_realnum, TScore_imgnum, TScore_phasenum, TScore_modnum, ref, nrows = 5)
    }
    return(result)
}
