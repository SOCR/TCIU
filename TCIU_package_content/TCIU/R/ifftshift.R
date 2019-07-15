#' @title ifftshift
#' @description IFFT SHIFT
#' @details This function is useful for moving back the zero-frequency component in the middle of the spectrum
#' back to (0,0,0).  It rearranges in reverse (relative to fftshift()) the indices appropriately,
#' so that the image can be correctly reconstructed by the IFT in spacetime
#'
#' @param img_ff An Inverse Fourier transform of a 1D signal, 2D image, or 3D volume.
#' @param dim Number of dimensions (-1, 1, 2, 3).
#' @return A properly shifted IFT of the input array.
#' @examples
#' ifftshift(fft(Real+1i*Imaginary, inverse = T))
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
ifftshift <- function(img_ff, dim = -1) {
  rows <- dim(img_ff)[1]
  cols <- dim(img_ff)[2]

  swap_up_down <- function(img_ff) {
    rows_half <- floor(rows/2)
    return(rbind(img_ff[((rows_half+1):rows), (1:cols)], img_ff[(1:rows_half), (1:cols)]))
  }
  swap_left_right <- function(img_ff) {
    cols_half <- floor(cols/2)
    return(cbind(img_ff[1:rows, ((cols_half+1):cols)], img_ff[1:rows, 1:cols_half]))
  }
  if (dim == -1) {
    img_ff <- swap_left_right(img_ff)
    return(swap_up_down(img_ff))
  }
  else if (dim == 1) {
    return(swap_up_down(img_ff))
  }
  else if (dim == 2) {
    return(swap_left_right(img_ff))
  }
  else if (dim == 3) {
    # Use the `abind` package to bind along any dimension a pair of multi-dimensional arrays
    # install.packages("abind")
    library(abind)
    planes <- dim(img_ff)[3]
    rows_half <- floor(rows/2)
    cols_half <- floor(cols/2)
    planes_half <- floor(planes/2)

    img_ff <- abind(img_ff[1:rows, 1:cols, ((planes_half+1):planes)],
                    img_ff[1:rows, 1:cols, 1:planes_half], along=3)
    img_ff <- abind(img_ff[1:rows, ((cols_half+1):cols), (1:planes)],
                    img_ff[1:rows, 1:cols_half, (1:planes)], along=2)
    img_ff <- abind(img_ff[((rows_half+1):rows), (1:cols), (1:planes)],
                    img_ff[(1:rows_half), (1:cols), (1:planes)], along=1)
    return(img_ff)
  }
  else {
    stop("Invalid dimension parameter")
  }
}
