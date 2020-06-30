# @title fmri_split_ab_bl
#
# @description A way to create continuous hexadecimal color vector for data that would be applied to
# \code{plotly} \strong{Notice:} This function mainly serves as a functional part for other functions. So we seldom apply this
# function by itself.  
# @param vect a vector or a list contains the data that would be referred to create the hexadecimal color vector 
# @param option control the function of this function. If choose 'auto', then this function will lead you to key in the
# space (x,y,z) parameters and time (time) parameter for this function to generate graphs.  
# 
# @details The function
# \code{fmri_split_ab_bl} is used to create continuous hexadecimal color vector for data.  
# 
# @author SOCR team <\url{http://socr.umich.edu/people/}> 
# 
# @return a vector or list representing the hexadecimal color
# @export 
# @examples

fmri_split_ab_bl <- function(vect, option = "vector") {
    if (option == "list") {
        overalllen <- length(which(vect != 0)) + 1
        ab_len <- length(which(vect > 0)) + 1
        bl_len <- length(which(vect < 0)) + 1
        nulllist <- as.list(1:overalllen)
        s <- seq_gradient_pal("#FFFF00", "#FFFFFF")(seq(0, 1, length.out = bl_len))
        for (i in 1:bl_len) {
            nulllist[[i]] <- c((i - 1)/overalllen, s[i])
        }
        s <- seq_gradient_pal("#FFFFFF", "#0000FF")(seq(0, 1, length.out = ab_len))
        for (i in 1:ab_len) {
            nulllist[[bl_len - 1 + i]] <- c((bl_len - 1 + i - 1)/overalllen, s[i])
        }
    } else if (option == "vector") {
        overalllen <- length(which(vect != 0)) + 1
        ab_len <- length(which(vect > 0)) + 1
        bl_len <- length(which(vect < 0)) + 1
        nulllist <- rep(NA, overalllen)
        s <- seq_gradient_pal("#FFFF00", "#FFFFFF")(seq(0, 1, length.out = bl_len))
        for (i in 1:bl_len) {
            nulllist[i] <- s[i]
        }
        s <- seq_gradient_pal("#FFFFFF", "#0000FF")(seq(0, 1, length.out = ab_len))
        for (i in 1:ab_len) {
            nulllist[[bl_len - 1 + i]] <- s[i]
        }
    }
    return(nulllist)
}
