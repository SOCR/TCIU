#' @title tensor-on-tensor regression on region of interest(ROI) of the brain
#' @description This function takes a 4d fMRI data and detects locations where stimulus occurs
#' on each region of interest(ROI) of the brain using \code{MultiwayRegression}. This function could be used as
#' an intermediate step of a three-phase analytics protocol to detect motor areas. The functions to implement this 
#' three-phase protocol in a consecutive order is \code{fmri_ROI_phase2}, \code{fmri_ROI_phase3} and \code{fmri_post_hoc} respectively.
#'
#' @param fmridata a 4d array which contains the spatial and temporal record of fmri result.
#' @param label_mask a 3d nifti or 3d array of data that shows the labeled brain atlas.
#' @param label_dict a dataframe or array or matrix to specify the indices and corresponding
#' names of the ROI. The input of this parameter could take one of the list outputs of the \code{fmri_ROI_phase2} function as a following step.
#' @param stimulus_idx a vector of the start time points of the time period when the fMRI data receives stimulation.
#' @param stimulus_dur a vector of the time period when the fMRI data receives stimulation.
#' @param fmri.design_order a parameter to specify the order of the polynomial drift terms in \code{fmri.design} function.
#' @param fmri.stimulus_TR a parameter to specify the time between scans in seconds in \code{fmri.stimulus} function.
#' @param rrr_rank a parameter to specify the assumed rank of the coefficient array in \code{rrr} function.
#' @param method a string that represents method for calculating p-values from tensor-on-tensor regression coefficients. 
#' There are 2 options: 't_test' and 'corrected_t_test'. The default is 't_test'.
#' 't_test' is to calculate the test statistics 't-value' across all voxels in the bounding box of ROI; 
#' 'corrected_t_test' is to calculate the test statistics 't-value' by first across each voxel on a temporal basis,
#' and then across all voxels in the bounding box of ROI.
#' @param parallel_computing a logical parameter to determine whether to use parallel computing to speed up the function or not.
#' The default is FALSE.
#' @param ncor number of cores for parallel computing. The default is the number of cores of the computer minus 2.
#' 
#' @details
#' The function \code{fmri_ROI_phase2} is used to detect locations where stimulus occurs by calculating the p-values
#' of the ROI-based tensor-on-tensor regression. Two methods can be chosen to calculate the p-values from the regression coefficients.
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return a 3d array storing ROI-based tensor regression p-values for the 4d fMRI data
#'
#' @examples
#' \donttest{
#' # sample 3D data of labeled brain atlas provided by the package
#' # this example will use parallel computing and take about ten minutes to finish
#' dim(mask_label)
#' # sample dataframe of ROI-based indices and names provided by the package
#' dim(mask_dict)
#' # sample 3D data of mask provided by the package
#' dim(mask)
#'
#' # calculated p-values 
#' set.seed(1)
#' fmri_generate = fmri_simulate_func(dim_data = c(64, 64, 40), mask = mask)
#' fmridata = fmri_generate$fmri_data
#' stimulus_idx = fmri_generate$ons
#' stimulus_dur = fmri_generate$dur
#' # the function will may take a long time, see examples in demo function or vignettes       
#' }
#'
#' @export
#'
#' @import dplyr fmri MultiwayRegression doParallel parallel
#' @importFrom stats pt sd


fmri_ROI_phase2 = 
  function(fmridata,
           label_mask,
           label_dict,  
           stimulus_idx,
           stimulus_dur,
           fmri.design_order = 2, 
           fmri.stimulus_TR = 3,
           rrr_rank = 3,
           method = "t_test",
           parallel_computing = FALSE, 
           ncor = max(detectCores()-2,1)
  ){

    # generate customized error message based on different wrong input
    if(length(dim(fmridata)) != 4){
      stop("'fmridata' should be a 4d array.")
    }
    
    if(((class(label_mask) %in% c("array", "nifti")) != TRUE) | (length(dim(label_mask)) != 3)){
      stop("'mask' should be a 3d array.")
      }
    else if(all.equal(dim(label_mask),dim(fmridata)[c(1,2,3)]) != TRUE){
      stop("The dimension of 'label_mask' should be the same as the first three dimension of 'fmridata'.")
    }
    
    if(!class(label_dict) %in% c("array","matrix","data.frame")){
      stop("'label_dict' should be an array or matrix or dataframe.")
    }
    else if(ncol(label_dict) < 2){
      stop("'label_dict' should have at least two columns as indices and names of the ROI.")
    }
    else if((class(label_dict[,1]) %in% c("integer","numeric") & class(label_dict[,2]) %in% c("character","factor")) != TRUE){
      stop("'label_dict' should have the first column as numeric indices and the second column as names of the ROI.")
    }
    
    # bounding box for ROI
    ROI_bounding_box = function(fmridata, label_mask, label_id){
      
      fmridata = Mod(fmridata) 
      time_span = dim(fmridata)[4]
      
      ROI_index = as.data.frame(which(label_mask == label_id, arr.ind = T))
      names(ROI_index) = c("x","y","z")
      x.min <- min(ROI_index$x)
      x.max <- max(ROI_index$x)
      y.min <- min(ROI_index$y)
      y.max <- max(ROI_index$y)
      z.min <- min(ROI_index$z)
      z.max <- max(ROI_index$z)
      
      ROI_index_move = as.data.frame(cbind(ROI_index[,1]-x.min+1,ROI_index[,2]-y.min+1,ROI_index[,3]-z.min+1))
      names(ROI_index_move) <- c("x","y","z")
      
      dim_for_block = c(x.max-x.min+1,y.max-y.min+1,z.max-z.min+1)
      bounding_box = array(0,dim = c(dim_for_block,time_span))
      for(t in 1:time_span){
        bounding_box[as.matrix(cbind(ROI_index_move,t))] = fmridata[as.matrix(cbind(ROI_index,t))]
      }
      
      list(bounding_box,ROI_index,ROI_index_move)
    }
    
    # p_value for ROI
    block_p_value = function(BOLD_coef, time_span, num_of_predictors){
      n = time_span
      dim_for_block = dim(BOLD_coef)
      p = num_of_predictors
      p_value = array(NA,dim_for_block)
      t_value = BOLD_coef/sd(BOLD_coef)
      p_value =  2 * pt(abs(t_value), df = n-p-1, lower.tail = FALSE)
      p_value
    }
    
    # main function
    time_span = dim(fmridata)[4]
    fixed_stim = fmri.stimulus(time_span, onsets = stimulus_idx, durations = stimulus_dur, TR = fmri.stimulus_TR)
    X_tensor = fmri.design(fixed_stim, order = fmri.design_order) 
    dim(X_tensor) = c(time_span, fmri.design_order + 2, 1)
    
    overall_p_value = array(1,dim(fmridata)[c(1,2,3)])
    label_list = label_dict[,1]
    label_name = as.character(label_dict[,2])
    
    switch(method, 
           "t_test" = {
             if(parallel_computing == TRUE){
               cl = makeCluster(ncor)
               registerDoParallel(cl)
               
               V = list()
               V = foreach(i = 1:length(label_list)) %dopar% {
                 label_id = label_list[i]
                 ROI =  ROI_bounding_box(fmridata, label_mask, label_id)
                 Y_tensor = aperm(ROI[[1]],c(4,1,2,3))
                 BOLD_coef = MultiwayRegression::rrr(X_tensor, Y_tensor, R = rrr_rank)$B[1,,,,]
                 p_value = block_p_value(BOLD_coef, time_span, fmri.design_order+2)
                 list(p_value,ROI)
               }
               stopCluster(cl)
               for(i in 1:length(label_list)){
                 overall_p_value[as.matrix(V[[i]][[2]][[2]])] = V[[i]][[1]][as.matrix(V[[i]][[2]][[3]])]
               }
             }
             
             
             else{
               for(i in 1:length(label_list)){
                 label_id = label_list[i]
                 ROI =  ROI_bounding_box(fmridata, label_mask, label_id)
                 Y_tensor = aperm(ROI[[1]],c(4,1,2,3))
                 BOLD_coef = MultiwayRegression::rrr(X_tensor, Y_tensor, R = rrr_rank)$B[1,,,,]
                 p_value = block_p_value(BOLD_coef, time_span, fmri.design_order+2)
                 overall_p_value[as.matrix(ROI[[2]])] = p_value[as.matrix(ROI[[3]])]
               }
             }
           },
           
           "corrected_t_test" = {
             if(parallel_computing == TRUE){
               cl = makeCluster(ncor)
               registerDoParallel(cl)
               
               V = list()
               V = foreach(i = 1:length(label_list)) %dopar% {
                 label_id = label_list[i]
                 ROI =  ROI_bounding_box(fmridata, label_mask, label_id)
                 Y_tensor = aperm(ROI[[1]],c(4,1,2,3))
                 BOLD_coef = MultiwayRegression::rrr(X_tensor, Y_tensor, R = rrr_rank)$B[1,,,,]
                 
                 # corrected p_value
                 sd_block = apply(ROI[[1]], c(1,2,3), sd)
                 corrected_BOLD = BOLD_coef/sd_block
                 corrected_BOLD[which(!is.finite(corrected_BOLD))] = 0
                 
                 p_value = block_p_value(corrected_BOLD, time_span, fmri.design_order+2)
                 list(p_value,ROI)
               }
               stopCluster(cl)
               
               for(i in 1:length(label_list)){
                 overall_p_value[as.matrix(V[[i]][[2]][[2]])] = V[[i]][[1]][as.matrix(V[[i]][[2]][[3]])]
               }
             }
             
             else{
               for(i in 1:length(label_list)){
                 
                 label_id = label_list[i]
                 ROI =  ROI_bounding_box(fmridata, label_mask, label_id)
                 Y_tensor = aperm(ROI[[1]],c(4,1,2,3))
                 BOLD_coef = MultiwayRegression::rrr(X_tensor, Y_tensor, R = rrr_rank)$B[1,,,,]
                 
                 # corrected p_value
                 sd_block = apply(ROI[[1]], c(1,2,3), sd)
                 corrected_BOLD = BOLD_coef/sd_block
                 corrected_BOLD[which(!is.finite(corrected_BOLD))] = 0
                 
                 p_value = block_p_value(corrected_BOLD, time_span, fmri.design_order+2)
                 overall_p_value[as.matrix(ROI[[2]])] = p_value[as.matrix(ROI[[3]])]

               }
             }
             
           }
     
    )
    
    overall_p_value
  }
