#' @title inverse kimesurface transform on a function in different periodic ranges
#' @description This function applies the inverse kimesurface transform to convert a kimesurface-transformed function back
#' to get the original 1D function in [0, 2*pi] or other similar periodic time range.
#' 
#' @param time_points a sequence of points in [0, 2*pi] or other periodic range
#' @param array_2d 2D array, got from the kimesurface_transform
#' @param num_length integer, interpolate f(t) to num_length samples in [0 : 2*pi] to extend the plot
#' @param m width of the contour path in C; too small values may lead to singularities on the negative x-axis; 
#' too large valued may lead to numerical instability for large positive x-axis. The default is 1.
#' @param msg Boolean to show/hide warnings. The default is TRUE.
#' 
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return a list of two elements
#' \itemize{
#'   \item Smooth_Reconstruction - the smoothed data computed from inverse kimesurface transform, 
#'   with the same length of time_points
#'   \item Raw_Reconstruction - the original unsmoothed data computed from inverse kimesurface transform, 
#'   with the same length of time_points
#' }
#' 
#' @examples
#' require(reshape2)
#' require(ggplot2)
#' \donttest{
#' # drop the first row and first column because of divergence on Laplace Transform
#' x = seq(0, 2, length.out=50)[2:50]; y = seq(0, 2, length.out=50)[2:50];
#' # do kimesurface transform on sine function
#' z2_grid = kimesurface_transform(FUNCT = function(t) { sin(t) },
#'                                 real_x = x, img_y = y)
#'                                 
#' time_points = seq(0+0.001, 2*pi, length.out = 160)
#' inv_data = inv_kimesurface_transform(time_points, z2_grid)
#' time_Intensities_ILT_df2 <- as.data.frame(cbind(Re=scale(Re(inv_data$Smooth_Reconstruction)), 
#'                                                 Im=scale(Re(inv_data$Raw_Reconstruction)),
#'                                                 fMRI=scale(Re(sin(time_points))),
#'                                                 time_points=time_points))
#' colnames(time_Intensities_ILT_df2) = c("Smooth Reconstruction", 
#'                                        "Raw Reconstruction", 
#'                                        "Original sin()",
#'                                        "time_points")
#' df = reshape2::melt(time_Intensities_ILT_df2, id.var = "time_points")
#' ggplot(df, aes(x = time_points, y = value, colour = variable)) + 
#'        geom_line(linetype=1, lwd=3) +
#'        ylab("Function Intensity") + xlab("Time") +
#'        theme(legend.position="top")+
#'        labs(title= bquote("Comparison between" ~ "f(t)=sin(t)" ~ "
#'        and Smooth(ILT(LT(fMRI)))(t); Range [" ~ 0 ~":"~ 2*pi~"]"))
#' }
#' @export
#' @importFrom reshape2 melt
#' @importFrom stats smooth.spline


inv_kimesurface_transform = function(time_points,
                                     array_2d,
                                     num_length = 20,
                                     m = 1,
                                     msg = TRUE){

  # this function convert the modified function value back to the original one
  inv_kimesurface_fun <- function (z, array_2D) {
    # convert z in C to Cartesian (x,y) coordinates
    x1 <- ceiling(Re(z))-1; # if (x1<2 || x1>dim(array_2d)[1]) x1 <- 2
    y1 <- ceiling(Im(z))-1; # if (y1<2 || y1>dim(array_2d)[2]) y1 <- 2
    # if exceed the domain use the default 1
    if(!is.na(x1)){
      if((x1 < 1) || (x1 > dim(array_2D)[1])){ x1 <- 1 }
    }
    if(!is.na(y1)){
      if((y1 < 1) || (y1 > dim(array_2D)[2])){ y1 <- 1 }
    }
    # exponentiate to Invert the prior (LT) log-transform of the kimesurface magnitude
    val1 = complex(length.out=1, real=Re(array_2D[x1, y1]), imaginary = Im(array_2D[x1, y1]))
    mag = exp(sqrt( Re(val1)^2+ Im(val1)^2))
    phase = atan2(Im(val1), Re(val1))
    value <- complex(real=Re(mag * exp(1i*phase)), imaginary = Im(mag * exp(1i*phase)))
    return ( value )
  }
  
  gamma=0.5
  fail_val = complex(1)
  nterms = 31L
  f2 <- array(complex(1), dim = length(time_points))
  for(t in 1:length(time_points)){
    
    f2[t] = ILT(FUNCT=function(z) inv_kimesurface_fun(z, array_2D=array_2d),
                 t= time_points[t], nterms,
                m, gamma, fail_val, msg)
  }

  # interpolate f(t) to 20 samples in [0 : 2*pi]. Note that 20~2*PI^2
  tvalsn <- seq(0, pi*2, length.out = num_length)
  f3 <- array(complex(1), dim = length(time_points)); # length(f), f=ILT(F)
  for (t in 1:length(time_points)) {
    f3[t] <- f2[ceiling(t/num_length)]
  }
  f4 <- smooth.spline(time_points, Re(f3), spar=1)$y; # plot(f4)
  
  return(list(Smooth_Reconstruction=f4, Raw_Reconstruction=f3))

}
