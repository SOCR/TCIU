#' @title p-values on region of interest(ROI) of the brain
#' @description This function takes a 4 dimensional real-valued fMRI data and calculates p-values for the ROIs individually to test whether the ROI is potentially activated.
#' It is the first phase of a ROI 3-phase analysis and usually followed by second phase analysis \code{fmri_ROI_phase2} .
#' @param fmridata a 4d array which contains the spatial and temporal record of fmri data
#' @param label_mask a 3D nifti or 3D array of data to indicates the corresponding indices of the ROIs
#' @param label_dict a dataframe which contains the name of ROIs and their corresponding index
#' @param stimulus_idx a vector that specifies when motion happens
#' @param rest_idx a vector that specifies when study participant does not move
#' @param p_threshold NULL or a numeric value that can be selected randomly below 0.05 to drop all p-values above the threshold.
#' @details
#' The function \code{fmri_ROI_phase1} is used to calculate p-values of ROIs for a given real-valued fmridata. It first takes in the fmridata and corresponding mask. 
#' For a fixed region, the function will first compute Temporal Contrast-to-noise Ratio (tCNR) for each voxel in that region, which is the mean of 80 paired differences in intensity for "on" and "off" states divided by its standard deviation. Second, it will conduct t-test on all tCNRs of a fixed region to see there are significant changes for the ROI during the on and off period. Finally, it will use bonferroni correction to control significant level and select the ROIs with p-values under the significant level to enter next phase analysis. 
#' 
#'@author SOCR team <\url{http://socr.umich.edu/people/}>
#
#' @return a list of two elements
#' \itemize{
#'   \item all_ROI - the test result for all ROIs
#'   \item sign_ROI - the test result for significant ROIs
#' }
#' 
#' @export
#'
#'
#' @examples 
#' fmri_generate = fmri_simulate_func(dim_data = c(64, 64, 40), mask = mask)
#' # p-values for phase 1
#' \donttest{
#' result = fmri_ROI_phase1(fmri_generate$fmri_data, mask_label,
#'                           mask_dict, stimulus_idx = fmri_generate$on_time)
#' }



fmri_ROI_phase1 = function(fmridata, 
                         label_mask=NULL,
                         label_dict=NULL,
                         stimulus_idx = NULL,
                         rest_idx = NULL,
                         p_threshold=0.05
                         
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
     stop("'ROI_label_dict' should be an array or matrix or dataframe.")
   }
   else if(ncol(label_dict) < 2){
     stop("'ROI_label_dict' should have at least two columns as indices and names of the ROI.")
   }
   else if((class(label_dict[,1]) %in% c("integer","numeric") & class(label_dict[,2]) %in% c("character","factor")) != TRUE){
     stop("'ROI_label_dict' should have the first column as numeric indices and the second column as names of the ROI.")
  }
  
  
  time_span = dim(fmridata)[4]
  
  
  ON_idx = stimulus_idx
 
  OFF_idx  = setdiff(1:time_span,ON_idx)
  
  
  if (is.null(rest_idx) == FALSE){
    OFF_idx = rest_idx
  }
  index=label_dict$index
  Y_ON_OFF=fmridata[,,,ON_idx]-fmridata[,,,OFF_idx]
 
  Y_ON_OFF_mat=matrix(Y_ON_OFF, prod(dim(Y_ON_OFF)[1:3]),dim(Y_ON_OFF)[4])
  mean_ON_OFF<- rowMeans(Y_ON_OFF_mat)
  sd_ON_OFF<-apply(Y_ON_OFF_mat,1,sd)

  
  CNR<-mean_ON_OFF/sd_ON_OFF
  CNR<-array(CNR,dim=dim(Y_ON_OFF)[1:3])
  pval_t<-vector(length=length(index))
  for (j in 1:length(index))
  {
    
    pval_t[j]=t.test(CNR[label_mask==index[j]])$p.value
  }
  all_ROI<-data.frame(index=index,name=label_dict$name, pval_t=pval_t)
  all_ROI<-all_ROI[order(all_ROI$pval_t),]
  sign_ROI<-all_ROI[all_ROI$pval_t<=p_threshold/length(index),]
  result<-list(all_ROI=all_ROI,sign_ROI=sign_ROI)
}

