# @title fmri_complex_p_val
# 
# @description This function takes a 4 dimensional complex-valued fmridata or a single complex valued
# vector, and calculates the p value for the given data to test whether there are significant changes for the given voxel during
# the on and off period.  
# 
# @param fmridata a 4d array which contains the spatial and temporal record of fmri result or a single
# complex-valued vector.  
# @param voxel_location a 3d array indicating the spatial location of the brain.  
# @param stimulus_idx a vector that specifies when motion happens.  
# @param rest_idx a vector that specifies the time points when study participant does not move.  
# @param method a string that indicates which testing method is to be used. There are 3 options: 'HotellingT2',
# 'Wilk's Lambda' and 'gLRT'(generalized likelihood ratio test). By default is 'HotellingT2'.  
# @param is.4d By default is true. If change to false, need to input a vector instead of a 4d array.  
# @param plot By default is false. If true, one can also see the visualization for the test.  
# 
# @details The function \code{complexfmri.p_val} is used to calculate p value for a given
# voxel of a complex-valued fmridata. It first takes in the fmridata, and then users need to specify whether the input is the
# whole array with the location for the voxel or the input is already a vector storing data for the given voxel. Users can then
# choose to use various methods to conduct hypothesis test to see whether there are significant changes for the given voxel
# during the on and off period. 
# 
# @author SOCR team <\url{http://socr.umich.edu/people/}> 
# 
# @return the p-value result for the given voxel of fmridata or the given voxel.  
# 
# @export 
# 
# @examples 
# # require('ICSNP') # require('rrcov') # simulate the fmridata
# fmri_generate = fmri_simulate_func(dim_data = c(64, 64, 40), mask = mask) p_val = complexfmri.p_val(fmri_generate$fmri_data,
# voxel_location = c(40,22,33), plot = TRUE)
fmri_complex_p_val = function(fmridata, voxel_location = NULL, method = "HotellingsT2", stimulus_idx = NULL, rest_idx = NULL, is.4d = TRUE, 
    plot = FALSE, ons = NULL, dur = NULL) {
    
    # Gives flexibility of either taking in a 4D array or a vector
    if (is.4d == TRUE && is.null(voxel_location) == FALSE) {
        voxel = fmridata[voxel_location[1], voxel_location[2], voxel_location[3], ]
    } else {
        voxel = fmridata
    }
    
    if (method == "gLRT") {
        return(complex_gen_est_pval(voxel = voxel, onsets = ons, durations = dur))
    }
    
    time_span = length(voxel)
    
    ON_idx = stimulus_idx
    # Unless specified, off_idx would be the complement of on_idx
    OFF_idx = setdiff(1:time_span, ON_idx)
    if (is.null(rest_idx) == FALSE) {
        OFF_idx = rest_idx
    }
    
    # Give labels to all the time series data
    labels = rep(1, time_span)
    labels[OFF_idx] = 0
    labels = factor(labels)
    
    Y1 = cbind(Re(voxel[ON_idx]), Im(voxel[ON_idx]))
    Y2 = cbind(Re(voxel[OFF_idx]), Im(voxel[OFF_idx]))
    # bind the 2 bivariate datasets into one DF
    Y12 <- rbind(Y1, Y2)  # View(cbind(Y12, labels))
    
    if (plot == TRUE) {
        # Visually inspect the 2D bivariate sim data
        plot(Y1[, 1], Y1[, 2], pch = "o", col = "blue", xlim = c(min(Y1[, 1], Y2[, 1]), max(Y1[, 1], Y2[, 1])), ylim = c(min(Y1[, 
            2], Y2[, 2]), max(Y1[, 2], Y2[, 2])), xlab = "dim1", ylab = "dim2", main = "Bivariate observations: timepoints=8*2*10=160 Hotelling's T2")
        graphics::points(Y2[, 1], Y2[, 2], pch = "*", col = "red")
    }
    switch(method, HotellingsT2 = {
        test = ICSNP::HotellingsT2(Y12 ~ labels)
    }, `Wilk's Lambda` = {
        test = rrcov::Wilks.test(Y12, grouping = labels, method = "c")
    })
    return(test)
}
