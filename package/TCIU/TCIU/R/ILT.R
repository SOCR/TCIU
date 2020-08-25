#' @title numerical method to compute inverse of Laplace Transform 
#' @description a function that numerically computes the inverse of Laplace Transform
#' 
#' @param FUNCT function object F(z), typically a Laplace Transform of a function f(t)
#' @param t time domain value to evaluate the ILT(F)(t)
#' @param nterms number of terms to use in the numerical inversion (odd number). The default is 31L.
#' @param m width of the contour path in C; too small values may lead to singularities on the negative x-axis; 
#' too large valued may lead to numerical instability for large positive x-axis. The default is 1.
#' @param gamma value on the positive x-axis for the vertical line representing the contour. The default is 0.5
#' @param fail_val value to return in event of failure to converge
#' @param msg Boolean to show/hide warnings. The default is TRUE.
#' 
#' @details
#' This function first uses full optimum contour path to do inverse Laplace Transform.
#' However, if this method fails, the function will automatically change to the method 
#' of using Bromwich contour path to do inverse Laplace Transform
#' 
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return a real value computed from inverse Laplace Transform
#' 
#' @examples
#' # analytic form of Laplace transform of f(t) = t
#' F = function(z) { 1/(z^2) }
#' # do inverse Laplace transform on t = 0.2
#' ILT(F, t = 0.2)
#' # the result is equal to t = 0.2
#' @export

ILT = function(FUNCT, t, nterms = 31L, m = 1, gamma = 0.5, fail_val = complex(0), msg = TRUE) {
    # some small built-in function for ILT 1 Cartesian to Polar Coordinate Transform
    r.xy <- function(x, y) {
        return(sqrt(x^2 + y^2))
    }
    
    phi.xy <- function(x, y) {
        return(atan2(y, x))
    }
    
    # 2 Polar to Cartesian Coordinate Transform
    x.rphi <- function(r, phi) {
        return(r * cos(phi))
    }
    
    y.rphi <- function(r, phi) {
        return(r * sin(phi))
    }
    
    # 3 Optimum contour path in Complex plane Uses Evans, Chung (2000) method
    optimContour <- function(phi, m = 1, t = 5) {
        if (identical(t, 0)) 
        {
        t <- 5
        }  # reset t=0 to avoid singularities
        return(m * phi/(t * sin(phi)))
    }
    
    analyticalPathDerivative <- function(phi, m = 1, t = 5) {
        if (identical(t, 0)) 
        {
        t <- 100
        }  # reset t = 0 to avoid singularities
        if (identical(phi, 0)) 
        {
        dx_dphi <- 0
        }  # avoid small phi value singularities
        else {
        dx_dphi <- (m/t) * ((sin(phi) - phi * cos(phi)) * cos(phi)/(sin(phi))^2 - phi)
        }
        dy_dphi <- (m/t)
        return(complex(real = dx_dphi, imaginary = dy_dphi))
    }
    
    # 4 bromwich contour path in complex plane Uses Christopher, Barry (2015) method
    bromwichContour = function(phi, gamma = 0.5) {
        return(gamma/cos(phi))
    }
    
    bromwichPathDerivative = function(phi, gamma = 0.5) {
        return((0 + (0 + (0 + (0+1i)))) * (bromwichContour(phi, gamma)^2 + (gamma * tan(phi)/cos(phi))^2)^0.5)
    }
    
    
    # main function
    FUNCT <- match.fun(FUNCT)  # input C-valued FUNCT should be interpreted as a function
    if (identical(t, 0)) 
    {
        t <- 10^-50
    }  # accurate for very low values of t>0
    n_attempts <- 10
    
    for (attempt in 1:n_attempts) {
        # repeated numerical integration attempts
        dphi <- 2 * pi/nterms
        phi <- -pi + (1:nterms - 1/2) * dphi
        z <- optimContour(phi, m, t) * cos(phi) + optimContour(phi, m, t) * sin(phi) * (0 + (0 + (0+1i)))
        # Plot optimal contour # plot(Re(z), Im(z))
        # Specify the return Function type must be a single C-value!
        # L <- sapply(z, FUNCT) # Funct returns array or list, not a single C-value!
        
        if (any(is.nan(z))) {
        #if original method produce nan, directly turn to 2nd method
        if (msg == TRUE & identical(attempt, as.integer(n_attempts))) {
            # if original method fails, construct ILT ob bromwich contour
            sprintf("Laplace Transform inversion failed after %s attempts.\n", n_attempts)
            # Note: Try alternative bromwich contour in C, which may not be very stable
            nterms = 1000L
            tot <- 0
            dphi <- pi/nterms
            for (n in 1:nterms) {
            phi <- -pi/2 + (n - 1/2) * dphi
            ds <- dphi * bromwichPathDerivative(phi, gamma)
            x <- x.rphi(bromwichContour(phi, gamma), phi)
            y <- y.rphi(bromwichContour(phi, gamma), phi)
            tot <- tot + exp((x + y * (0 + (0 + (0 + (0+1i))))) * t) * FUNCT(x + y * (0 + (0 + (0 + (0+1i))))) * ds
            }
            return(tot/(2 * pi * (0 + (0 + (0 + (0+1i))))))
        }
        m <- m * 2  # contour path increment
        nterms <- round(nterms * 1.378, 0)  # random irregular number to avoid sampling same point again
        next
        }
        else
        {
        L <- vapply(z, FUNCT, 0 + (0 + (0+1i))) 
        if (any(!is.finite(L))) {
            if (msg == TRUE & identical(attempt, as.integer(n_attempts))) {
            # if original method fails, construct ILT ob bromwich contour
            sprintf("Laplace Transform inversion failed after %s attempts.\n", n_attempts)
            # Note: Try alternative bromwich contour in C, which may not be very stable
            nterms = 1000L
            tot <- 0
            dphi <- pi/nterms
            for (n in 1:nterms) {
                phi <- -pi/2 + (n - 1/2) * dphi
                ds <- dphi * bromwichPathDerivative(phi, gamma)
                x <- x.rphi(bromwichContour(phi, gamma), phi)
                y <- y.rphi(bromwichContour(phi, gamma), phi)
                tot <- tot + exp((x + y * (0 + (0 + (0 + (0+1i))))) * t) * FUNCT(x + y * (0 + (0 + (0 + (0+1i))))) * ds
            }
            return(tot/(2 * pi * (0 + (0 + (0 + (0+1i))))))
            }
            m <- m * 2  # contour path increment
            nterms <- round(nterms * 1.378, 0)  # random irregular number to avoid sampling same point again
            next
        }
        }
        em <- exp(z * t)
        ds_dphi <- vapply(phi, analyticalPathDerivative, 0 + (0 + (0+1i)), m, t)
        break
    }
    return(sum(em * L * ds_dphi) * dphi/((0 + (0 + (0+2i))) * pi))
}
