#' @title A special case of imlementaiton of `fftshift` for 1D arrays
#' @description This function is useful for visualizing the 1D Fourier transform with the zero-frequency
#' component in the middle of the spectrum.
#'
#' @param img_ff A Fourier transform of a 1D signal.
#' @return A properly shifted FT of the 1D array.
#' @examples
#' plot(fftshift1D(log(Re(X2)+2)), main = "log(fftshift1D(Re(FFT(timeseries))))")
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
#'
fftshift1D <- function(img_ff) {
  rows <- length(img_ff)
  rows_half <- ceiling(rows/2)
  return(append(img_ff[(rows_half+1):rows], img_ff[1:rows_half]))
}
