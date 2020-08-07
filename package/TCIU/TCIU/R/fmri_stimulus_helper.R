# Helper funtions used in fmri_stimulus_detect to calculate p-values

# Calculate p-value for a given voxel of a real valued fmridata
# Two options: "t-test" and "wilcoxon-test"
fmri_p_val = function(fmridata, 
                      voxel_location = NULL,
                      stimulus_idx = NULL,
                      rest_idx = NULL,
                      epoch_length = 10,
                      is.4d = TRUE,
                      test_type = "t-test"){

  # Get the averages for the 8 finger-tapping epochs
  # Give flexibility of either taking in a 4D array or a vector
  if (is.4d == TRUE & is.null(voxel_location) == FALSE ){
    voxel = fmridata[voxel_location[1],
                     voxel_location[2],
                     voxel_location[3], ]
  }else{
    voxel = fmridata
  }
  
  time_span = length(voxel)
  
  # Take into account the hemodynamic response delay
  # data_seq = c(data_seq, rep( mean( data_seq[153:160] ), 2) )
  
  ON_idx = stimulus_idx
  # Unless specified, off_idx would be the complement of on_idx
  OFF_idx  = setdiff(1:time_span,ON_idx)
  
  
  if (is.null(rest_idx) == FALSE){
    OFF_idx = rest_idx
  }
  
  # Get the averages for the 8 active-state epochs
  group1_avg = 
    voxel[ON_idx] %>%
    matrix(nrow = epoch_length) %>% 
    apply(2,mean)
  
  # Get the averages for the 8 rest-state epochs
  group2_avg = 
    voxel[OFF_idx] %>%
    matrix(nrow = epoch_length) %>%
    apply(2,mean)
  
  switch(test_type, 
         "t-test"= {
           return(t.test(group1_avg,group2_avg, 
                         alternative="greater")$p.value)},
         "wilcoxon-test"= {
           return( 
             wilcox.test(group1_avg, 
                         group2_avg,
                         alternative="greater" )$p.value )},
         
         {return("Please type a valid test type!!!")}
  )
}

# Calculate p-value for a given voxel of a real valued fmridata
# Three options: 'HotellingT2', 'Wilks-Lambda' and 'gLRT'
fmri_complex_p_val = function(fmridata,
                              voxel_location = NULL,
                              method = "HotellingsT2",
                              stimulus_idx = NULL,
                              rest_idx = NULL,
                              is.4d = TRUE,
                              ons = NULL,
                              dur = NULL) {
    
    # Give flexibility of either taking in a 4D array or a vector
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
    # Bind the 2 bivariate datasets into one DF
    Y12 <- rbind(Y1, Y2)  # View(cbind(Y12, labels))
    
    switch(method, 
          "HotellingsT2" = {
                            test = ICSNP::HotellingsT2(Y12 ~ labels)
                           }, 
          "Wilks-Lambda" = {
                            test = rrcov::Wilks.test(Y12, grouping = labels, method = "c")
                           }
          )
    return(test)
}

# Localize the position of motor area in human brain with linear regression against design matrix and then draw inference on corresponding coefficients
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


# Calculate volume of on-off-difference period polar.volume of a voxel.
fmri_on_off_volume = function(data,
                              x,
                              y,
                              z,
                              coordinates = "polar"){
  voxel_data = data[x,y,z,]
  On_data = voxel_data[rep(c(TRUE,FALSE), c(10,10)) ]
  Off_data = voxel_data[rep(c(FALSE,TRUE), c(10,10)) ]
  
  switch(coordinates, 
         "polar"= {
           volume.diff = 1/80 * ( sum( On_data^2 ) - sum( Off_data^2 ) )
           return(volume.diff)},
         
         "cartesian"= {
           # randomly generate 8 phi kime-phases for each of the 10 time radii
           phi_8_vec <- matrix(NA, ncol=10, nrow = 8)
           for (t in 1:10) { 
             # for a given t, generate 8 new phases
             set.seed(t);
             phi_8_vec[ ,t] <-
               extraDistr::rlaplace(8,mu=0,sigma=0.5)
             # rank-order the phases for consistency
             # within the same foliation leaf
             phi_8_vec[ ,t] <- sort(phi_8_vec[ ,t])
             # force phases in [-pi: pi)
             for (i in 1:8) {
               if (phi_8_vec[i,t] < -pi) 
                 phi_8_vec[i,t] <- -pi
               if (phi_8_vec[i,t] >= pi) 
                 phi_8_vec[i,t] <- pi
             }
           }
           matrix_ON <- matrix(0, nrow = 21, ncol = 21) 
           matrix_OFF <- matrix(0, nrow = 21, ncol = 21) 
           for (t in 1:10) {
             for (p in 1:8) {
               x = 11+t*cos(phi_8_vec[p,t])
               y = 11+t*sin(phi_8_vec[p,t])
               matrix_ON[x,y]  <- On_data[(p-1)*10 +t]
               matrix_OFF[x,y] <- Off_data[(p-1)*10 +t]
             }
           }
           volume.diff = 1/(21^2) * (sum( matrix_ON^2 ) - sum( matrix_OFF^2 ))
           return( volume.diff )},
         
         {return("Please type a valid coordinate system type!!!")}
  )
  
}
