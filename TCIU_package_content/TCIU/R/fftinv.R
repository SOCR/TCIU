#' @title fftinv
#' @description Implicitly Invert the FT (IFT)
#' @details This function does the IFT and scales appropriately the  result to ensure IFT(FT()) = I()
#'
#' @param x FT of a dataset.
#' @return The IFT of the input array.
#' @examples
#' fftinv(fft(circle_arr))
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
#'
fftinv <- function( x ){
  fft( x, inverse=TRUE ) / length( x )
}
