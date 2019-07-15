#' @title fftshift
#' @description FFT SHIFT
#' @details This function is useful for visualizing the Fourier transform with the zero-frequency
#' component in the middle of the spectrum.
#'
#' @param img_ff A Fourier transform of a 1D signal, 2D image, or 3D volume.
#' @param dim Number of dimensions (-1, 1, 2, 3).
#' @return A properly shifted FT of the array.
#'
#' @examples
#' fftshift(ft_circle)
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
#'
fftshift <- function(img_ff, dim = -1) {
  rows <- dim(img_ff)[1]
  cols <- dim(img_ff)[2]
  # planes <- dim(img_ff)[3]
  swap_up_down <- function(img_ff) {
    rows_half <- ceiling(rows/2)
    return(rbind(img_ff[((rows_half+1):rows), (1:cols)], img_ff[(1:rows_half), (1:cols)]))
  }
  swap_left_right <- function(img_ff) {
    cols_half <- ceiling(cols/2)
    return(cbind(img_ff[1:rows, ((cols_half+1):cols)], img_ff[1:rows, 1:cols_half]))
  }
  #swap_side2side <- function(img_ff) {
  #  planes_half <- ceiling(planes/2)
  #  return(cbind(img_ff[1:rows, 1:cols, ((planes_half+1):planes)], img_ff[1:rows, 1:cols, 1:planes_half]))
  #}
  if (dim == -1) {
    img_ff <- swap_up_down(img_ff)
    return(swap_left_right(img_ff))
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
    rows_half <- ceiling(rows/2)
    cols_half <- ceiling(cols/2)
    planes_half <- ceiling(planes/2)

    img_ff <- abind(img_ff[((rows_half+1):rows), (1:cols), (1:planes)],
                    img_ff[(1:rows_half), (1:cols), (1:planes)], along=1)
    img_ff <- abind(img_ff[1:rows, ((cols_half+1):cols), (1:planes)],
                    img_ff[1:rows, 1:cols_half, (1:planes)], along=2)
    img_ff <- abind(img_ff[1:rows, 1:cols, ((planes_half+1):planes)],
                    img_ff[1:rows, 1:cols, 1:planes_half], along=3)
    return(img_ff)
  }
  else {
    stop("Invalid dimension parameter")
  }
}
