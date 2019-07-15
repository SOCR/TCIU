#' @title K-space Transformation
#' @description This function is useful for transforming data to k-space.
#' @details A generic function to Transform Data ={all predictors (X) and outcome (Y)} to k-space (Fourier domain)
#' For ForwardFT, set parameters as (rawData, FALSE, NULL)
#' For InverseFT, there are two parameters setting: (magnitudes, TRUE, reconPhasesToUse) or (FT_data, TRUE, NULL)
#'
#' @param data dataset that needs K-space transformation
#'
#' @examples
#' kSpaceTransform(mag_ft_xC_fMRI_train, TRUE, phase_ft_xC_fMRI_train[1:length(xC_fMRI_train)])
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
#'
kSpaceTransform <- function(data, inverse = FALSE, reconPhases = NULL) {
  # ForwardFT (rawData, FALSE, NULL)
  # InverseFT(magnitudes, TRUE, reconPhasesToUse) or InverseFT(FT_data, TRUE, NULL)
  FT_data <- array(complex(), c(dim(data)[1], dim(data)[2]))
  mag_FT_data <- array(complex(), c(dim(data)[1], dim(data)[2]))
  phase_FT_data <- array(complex(), c(dim(data)[1], dim(data)[2]))
  IFT_reconPhases_data <- array(complex(), c(dim(data)[1], dim(data)[2]))

  for (i in 1:dim(data)[2]) {
    if (inverse == FALSE | is.null(reconPhases)) {
      FT_data[ , i] <- fft(data[ , i], inverse)
      X2 <- FT_data[ , i]
      # plot(fftshift1D(log(Re(X2)+2)), main = "log(fftshift1D(Re(FFT(timeseries))))")
      mag_FT_data[ , i] <- sqrt(Re(X2)^2+Im(X2)^2);
      # plot(log(fftshift1D(Re(mag_FT_MCSI_data))), main = "log(Magnitude(FFT(timeseries)))")
      phase_FT_data[ , i] <- atan2(Im(X2), Re(X2));
      # plot(Re(fftshift1D(phase_FT_MCSI_data[ , 1])), main = "Shift(Phase(FFT(timeseries)))")
    }
    else {  # for IFT synthesis using user-provided Phases, typically from kime-phase aggregators
      Real <- data[ , i] * cos(reconPhases[ , i])
      Imaginary <- data[ , i] * sin(reconPhases[ , i])
      IFT_reconPhases_data[ ,i] <-
        Re(fft(Real+1i*Imaginary, inverse = TRUE)/length(data[ , i]))
    }
  }
  ######### Test the FT-IFT analysis-synthesis back-and-forth transform process
  #         to confirm calculations
  # X2 <- FT_data[ , 1]; mag_FT_data[ , 1] <- sqrt(Re(X2)^2+Im(X2)^2);
  # phase_FT_data[ , 1] <- atan2(Im(X2), Re(X2));
  # Real2 = mag_FT_data[ , 1] * cos(phase_FT_data[ , 1])
  # Imaginary2 = mag_FT_data[ , 1] * sin(phase_FT_data[ , 1])
  # man_hat_X2 = Re(fft(Real2 + 1i*Imaginary2, inverse = T)/length(X2))
  # ifelse(abs(man_hat_X2[5] - data[5, 1]) < 0.001, "Perfect Syntesis", "Problems!!!")
  #########

  if (inverse == FALSE | is.null(reconPhases)) {
    return(list("magnitudes"=mag_FT_data, "phases"=phase_FT_data))
    # Use kSpaceTransform$magnitudes & kSpaceTransform$phases to retrieve teh Mags and Phases
  }
  else {
    return(IFT_reconPhases_data)
    # Use Re(kSpaceTransform) to extract spacetime Real-valued reconstructed data
  }
}
