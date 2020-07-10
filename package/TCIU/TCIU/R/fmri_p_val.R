# @title fmri_p_val
# @description This function takes a 4 dimensional real valued fmridata or a single real valued vector, calculates the p value for the given data to test whether there is significant changes for the given voxel during the on and off period.
#
# @param fmridata a 4d array which contains the spatial and temporal record of fmri result or a single real valued vector.
# @param voxel_location a 3d array indicating the spacial location of the brain.
# @param stimulus_idx a vector that specifies when motion happens.
# @param rest_idx a vector that specifies when study participant does not move.
# @param test_type a string that indicates which testing method is to be used. There are 2 options: "t-test" and "wilcoxon-test". By default is "t-test".
# @param is.4d By default is true. If change to false, need to input a vector instead of array.
#
# @details
# The function \code{fmri_p_val} is used to calculate p value for a given voxel of a real valued fmridata. It first takes in the fmridata, and then users need to specify whether the input is the whole array with the location for the voxel or the input is already a vector storing data for the given voxel. Users can then choose to use various methods to conduct hypothesis test to see whether there is significant changes for the given voxel during the on and off period.
#
# @author SOCR team <\url{http://socr.umich.edu/people/}>
#
# @return the test result for the given voxel of fmridata.
# @export
#
# @examples
# # simulate the fmridata
# fmri_generate = fmri_simulate_func(dim_data = c(64, 64, 40), mask = mask)
# p_val = fmri.p_val(fmri_generate$fmri_data, voxel_location = c(40,41,33))
# 
# 
fmri_p_val = function(fmridata, 
                      voxel_location = NULL,
                      stimulus_idx = NULL,
                      rest_idx = NULL,
                      epoch_length = 10,
                      is.4d = TRUE,
                      test_type = "t-test"){
  # Get the averages for the 8 finger-tapping epochs
  # data_seq = data_4d[x,y,z,]
  
  # Gives flexibility of either taking in a 4D array or a vector
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