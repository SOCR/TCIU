# @title fmri_hrf_p_val
#
# @description an fMRI motor area localization function, used to analyze fMRI data with linear regression
# against design matrix. 
#
# @param fmridata a 4D real-valued array that corresponds to the functional MRI scan result
# @param mask a 3D array data that is the mask of the MRI data. The mask is like the brain shell, which is the boundary of the MRI data. The
# default is NULL.  
# @param ons a vector with the first time points of the time periods when the fMRI data receives stimulation
# @param dur a vector of same length as 'ons'. The element inside represents the time length of each stimulated period. Notice
# that the unsimulated period has the same length of time period as the stimulated period happened just before it.  
#
# @details The function hrf_p_val is used to help localize the position of motor area in human brain.  with linear regression against design
# matrix and then draw inference on corresponding coefficients.  The p-values of the significance of coefficients would be
# stored as a 3D array. 
#
# @author SOCR team <\url{http://socr.umich.edu/people/}> 
# @return a 3D array of p-values 
#
# @export 
# @import stats dplyr fmri
fmri_hrf_p_val = function(fmridata, mask = NULL, ons, dur) {
    # Get the dimensions of the fMRI data and take the Mod
    dim1 = dim(fmridata)[1]
    dim2 = dim(fmridata)[2]
    dim3 = dim(fmridata)[3]
    dim4 = dim(fmridata)[4]
    fmridata_mod <- Mod(fmridata)
    # Make the design matrix in which 1st column is the study design, 2nd col is the linear term, 3rd col is the quadratic term...
    fixed_stim <- fmri.stimulus(dim4, onsets = ons, durations = dur, TR = 3)
    design_matrix <- fmri.design(fixed_stim)
    p_val_3d = array(1, dim = c(dim1, dim2, dim3))
    
    for (i in 1:dim1) {
        for (j in 1:dim2) {
            for (k in 1:dim3) {
                tryCatch({
                  p = 1
                  # If use inputs a mask
                  if (is.null(mask) == FALSE) {
                    if (mask[i, j, k] == 1) {
                      # Fit a linear regression against the design matrix where the reponse is the mod data
                      lm_data <- cbind(fmridata_mod[i, j, k, ], design_matrix) %>% as.data.frame()
                      colnames(lm_data) <- c("Y", "X1", "X2", "X3", "X4")
                      lm_result = summary(lm(Y ~ ., data = lm_data))
                      # Analyze the significance of the coefficients
                      p = lm_result$coefficients %>% as.data.frame() %>% .["X1", "Pr(>|t|)"]
                    }
                    # If no mask is used
                  } else {
                    # Fit a linear regression against the design matrix where the reponse is the mod data
                    lm_data <- cbind(fmridata_mod[i, j, k, ], design_matrix) %>% as.data.frame()
                    colnames(lm_data) <- c("Y", "X1", "X2", "X3", "X4")
                    lm_result = summary(lm(Y ~ ., data = lm_data))
                    # Analyze the significance of the coefficients
                    p = lm_result$coefficients %>% as.data.frame() %>% .["X1", "Pr(>|t|)"]
                  }
                  # p = as.numeric( (ICSNP::HotellingsT2(Y12 ~ labels))$p.value )
                }, error = function(e) {
                })  # end of try-catch
                p_val_3d[i, j, k] = p
            }
        }
        p_val_3d[i, , ][is.na(p_val_3d[i, , ])] = 1
    }
    
    return(p_val_3d)
}
