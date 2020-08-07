#' @title numerical method to compute Laplace Transform 
#' @description a function that numerically computes the Laplace Transform
#' 
#' @param FUNCT a function object f(t) conducting Laplace Transform
#' @param z a complex domain value used to evaluate the F(z)=LT(f)(z)
#' 
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return a complex value computed from Laplace Transform
#' 
#' @examples
#' f = function(t) { t }; z= 1+1i; 
#' LT(f, z); 
#' # compare with the result from analytic form of Laplace Transform of f(t) = t
#' # analytic form is below
#' F = function (z) { 1/z^2 }; F(z)
#' # the two results are the same
#' @export
#'
#' @import cubature

LT = function(FUNCT, z) {
    FUNCT <- match.fun(FUNCT)  
    # input R-domain FUNCT should be interpreted as a function
    # define the integrand function (complex-valued) Note that cubic numeric integration is used (Unified Cubature Integration
    # Interface), Integral limits are exact [0, Inf)
    integrand <- function(t) {
        return(FUNCT(t) * exp(-z * t))
    }
    return(cubintegrate(f = function(t) Re(integrand(t)), lower = 0, upper = Inf, method = "pcubature")$integral + (0 + (0 + (0+1i))) * 
        cubintegrate(f = function(t) Im(integrand(t)), lower = 0, upper = Inf, method = "pcubature")$integral)
}
