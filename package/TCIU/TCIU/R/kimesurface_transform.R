#' @title kimesurface transform on a function with a specified set of complex values
#' @description a function applies the kimesurface transform on a function with a specified set of complex values
#' 
#' @param FUNCT function object f(t) to conduct kimesurface transform on
#' @param glb_para a vector of global objections that needed to be imported when using parallel computing
#' @param real_x a list of numeric values, which is the real part of a set of complex values
#' @param img_y a list of numeric values, which is the imaginary part of the set of complex values stated above
#' @param parallel_computing logical object to determine whether to use parallel computing to speed up the function or not.
#' The default is FALSE.
#' @param ncor number of cores for parallel computing. The default is 6.
#' 
#' @details
#' This function applies the kimesurface transform on a 1D function f(t), to have it converted to a 2D function. The input
#' is a set of complex values with the same number of real and imaginary parts. These two parts can specify a 2D plane 
#' of the same length and width. The new 2D function is defined on this 2D plane. It mainly does a 
#' Laplace Transform and modifies all the function values in a specific way to have them looks better in the plot. 
#' 
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return a 2d array that did kimesurface transform for the set of complex value (the real and imaginary parts can
#' construct a 2d plane)
#' 
#' @examples
#' # drop the first row and first column because of divergence on Laplace Transform
#' # do kimesurface transform on sine function
#' x = seq(0, 2, length.out=50)[2:50]; y = seq(0, 2, length.out=50)[2:50];
#' 
#' \donttest{
#' kimesurface_transform(FUNCT = function(t) {sin(t)}, real_x = x, img_y = y);
#' }
#'                       
#' @export
#' 
#' @import doParallel foreach


kimesurface_transform = function(FUNCT,
                                 glb_para,
                                 real_x,
                                 img_y,
                                 parallel_computing = FALSE,
                                 ncor = 6){
  
  if(parallel_computing == TRUE){
    # kime surface transform
    # use parallel computing to speed up code
    glb_parals = c(glb_para, 'cubintegrate', 'LT')
    
    cl <- makeCluster(ncor)
    registerDoParallel(cl)
    F = list()
    for (i in 1:length(real_x) ){
      F[[i]] = 
        foreach(j = 1:length(img_y),
                .export=glb_parals, #'cubintegrate', # global object that need to be import in parallel computing
                .packages='cubature') %dopar% {
                  F_result = LT(FUNCT, complex(real=real_x[i], imaginary = img_y[j]))
                  mag = log(sqrt( Re(F_result)^2+ Im(F_result)^2))   
                  # log-transform the magnitude to temper the kimesurface amplitude
                  phase = atan2(Im(F_result), Re(F_result))
                  mag * exp(1i*phase)
                }
      }
    stopCluster(cl)
    F_vec = lapply(F, unlist)
    array_2d = unlist(do.call(rbind, F_vec))
    
  }else{
    F_result <- array(complex(1), dim = c(length(real_x), length(img_y)))
    for (i in 1:length(real_x) ){
        for(j in 1:length(img_y)){
          F_result[i, j] = LT(FUNCT, complex(real=real_x[i], imaginary = img_y[j]))
        }
      }
    mag = log(sqrt( Re(F_result)^2+ Im(F_result)^2))   
    # log-transform the magnitude to temper the kimesurface amplitude
    phase = atan2(Im(F_result), Re(F_result))
    array_2d = mag * exp(1i*phase)
  }
  
  return(array_2d)
  
}

