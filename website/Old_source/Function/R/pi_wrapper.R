#' @title pi_wrapper
#' @description Define a pi-wrapper calculator to ensure phase-arithmetic stays within [-pi : pi)
#' @details This function makes sure that phase-arithmetic is between [-pi, pi]
#'
#' @param x A numeric value
#' @return A number that satifies one of the 4 cases specified below
#' @examples
#' pi_wrapper(pi/8 + 1.186667 )
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
pi_wrapper <- function (x){
  if(abs(x%%(2*pi)) <= pi) return(x%%(2*pi))
  else if (-2*pi <= x%%(2*pi) & x%%(2*pi) < -pi) return (2*pi + x%%(2*pi))
  else if (pi <= x%%(2*pi) & x%%(2*pi) < 2*pi) return (x%%(2*pi) - 2*pi)
  else return (0)
}
