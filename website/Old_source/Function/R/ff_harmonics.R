#' @title ff_harmonics
#' @description Compute the spectral decomposition of an array (harmonics)
#' @details This function computes the FT of the singlal and plots the first few harmonics
#'
#' @param x Original signal (1D, 2D, or 3D array).
#' @param n Number of first harmonics to report (integer).
#' @param up Upsamping rate (default=10).
#' @param plot Boolean indicating whether to print the harmonics plot(default==TRUE).
#' @param add whether to overplot the harmonics on an existing graph (default=FALSE),
#' @param main Title for the plot.
#' @return A plot and a dataframe with the sampled harmonics and their corresponding FT magnitudes/amplitudes.
#' @examples
#' ff_harmonics(x = y, n = 12L, up = 100L, col = 2L, lwd=3, cex=2)
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
#'
ff_harmonics = function(x=NULL, n=NULL, up=10L, plot=TRUE, add=F, main=NULL, ...) {
  # The discrete Fourier transformation
  dff = fft(x)
  # time
  t = seq(from = 1, to = length(x))

  # Upsampled time
  nt = seq(from = 1, to = length(x)+1-1/up, by = 1/up)

  #New spectrum
  ndff = array(data = 0, dim = c(length(nt), 1L))
  ndff[1] = dff[1] # mean, DC component
  if(n != 0){
    ndff[2:(n+1)] = dff[2:(n+1)] # positive frequencies come first
    ndff[length(ndff):(length(ndff) - n + 1)] = dff[length(x):(length(x) - n + 1)] # negative frequencies
  }

  # Invert the FT
  indff = fft(ndff/length(y), inverse = TRUE)
  idff = fft(dff/length(y), inverse = TRUE)
  if(plot){
    if(!add){
      plot(x = t, y = x, pch = 16L, xlab = "Time", ylab = "Measurement",
           col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
           main = ifelse(is.null(main), paste(n, "harmonics"), main))
      lines(y = Mod(idff), x = t, col = adjustcolor(1L, alpha = 0.5))
    }
    lines(y = Mod(indff), x = nt, ...)
  }
  ret = data.frame(time = nt, y = Mod(indff))
  return(ret)
}
