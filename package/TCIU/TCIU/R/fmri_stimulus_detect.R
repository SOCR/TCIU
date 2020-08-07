#' @title fMRI data stimulus detection
#' @description This function takes a real/complex valued fMRI data and detects locations where stimulus occurs
#'
#' @param fmridata an array or a vector which contains the spatial and/or temporal record of fMRI result
#' @param mask a 3d array indicating the spatial location of the brain
#' @param stimulus_idx a vector that specifies when motion happens 
#' @param rest_idex a vector that specifies when study participant does not move
#' @param method a string that indicates which testing method is to be used. There are 5 options: 'HotellingT2', 'Wilks-Lambda' and 'gLRT'(likelihood ratio test) for complex fMRI data and 't-test', 'wilcoxon-test' for real fMRI data. For 4D real-valued fMRI data, two more options: 'on_off_diff' and 'HRF' method can be applied.
#' @param fdr_corr a logical variable. True if FDR correction is to be applied
#' @param spatial_cluster.thr threshold p-value to be used for spatial clustering
#' @param spatial_cluster.size number of spatially connected voxels to be tested for spatial clustering
#' @param ons a vector with the first time points of the time periods when the fMRI data receives stimulation. The default is NULL. Need to specify when choose the method 'gLRT' or 'HRF'.
#' @param dur a vector of the time length of each stimulated period. The default is NULL. Need to specify when choose the method 'gLRT' or 'HRF'.
#'
#' @details
#' The function \code{fmri_stimulus_detect} is used to conduct motor area detection. It first takes in a real or complex valued fMRI data, and then users can choose to use various methods to find the spatial regions where motor area is located inside the brain. User can either input the 4d fMRI data and get a 3d array storing p-values or input the fMRI data with smaller dimension (e.g. fix the x,y axis) and get a vector storing p-values. Besides, one can use this function to just calculate raw p-values, and we also provide options so that users can do FDR correction and spatial clustering to get a more accurate result.
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return If input fMRI data is 4d, return a 3d array storing p-values for the 4d fMRI data. If input fMRI data is less than 4d, return a vector storing p-values for the fMRI data.
#' @export
#'
#' @import AnalyzeFMRI dplyr fmri
#' @importFrom ICSNP HotellingsT2
#' @importFrom rrcov Wilks.test
#' @importFrom extraDistr rlaplace
#' @importFrom stats arima.sim ecdf ksmooth lm na.omit p.adjust pchisq qchisq quantile runif smooth t.test ts wilcox.test
#' @useDynLib TCIU, .registration = TRUE
#'
#' @examples
#' fmri_generate = fmri_simulate_func(dim_data = c(64, 64, 40), mask = mask)
#' fmridata = fmri_generate$fmri_data
#' stimulus_idx = fmri_generate$on_time
#' ons = fmri_generate$ons
#' dur = fmri_generate$dur
#' 
#' \donttest{
#' # p-values using t-test for 4d fMRI data
#' p_value1 = fmri_stimulus_detect(fmridata = fmridata, mask = mask,
#'                                 stimulus_idx = stimulus_idx,
#'                                 method = 't-test')
#' dim(fmridata)
#' dim(p_value1)
#' 
#' # p-values using t-test for 2d fMRI data
#' p_value2 = fmri_stimulus_detect(fmridata = fmridata[40,41,,], mask = mask,
#'                                 stimulus_idx = stimulus_idx,
#'                                 method = 't-test')
#' dim(fmridata[40,41,,])
#' dim(p_value2)
#' }
#' 
fmri_stimulus_detect = function(fmridata,
                                mask = NULL,
                                stimulus_idx = NULL,
                                rest_idex = NULL,
                                method,
                                fdr_corr = NULL,
                                spatial_cluster.thr = NULL,
                                spatial_cluster.size = NULL,
                                ons = NULL,
                                dur = NULL) {
    
    # if data is less than 4d
    if (is.null(dim(fmridata)) == TRUE) {
        # 1d
        if (typeof(fmridata) != "complex" & method %in% c("t-test", "wilcoxon-test")) {
            p_val = fmri_p_val(fmridata = fmridata, stimulus_idx = stimulus_idx, test_type = method)
            if (is.na(p_val) == TRUE) {
                p_val = 1
            }
            return(p_val)
        } else if (typeof(fmridata) == "complex" & method %in% c("HotellingsT2", "Wilks-Lambda", "gLRT")) {
            p_val = tryCatch({
                as.numeric(fmri_complex_p_val(fmridata = fmridata, stimulus_idx = stimulus_idx, method = method, ons = ons, dur = dur)$p.value)
            }, error = function(e) {
                return(1)
            })  # end of try-catch
            return(p_val)
        } else {
            # wrong test method input
            if (typeof(fmridata) == "complex") {
                stop("Invalid test type! You can only use HotellingsT2, Wilks-Lambda or gLRT for complex data!")
            } else {
                stop("Invalid test type! You can only use t-test or wilcoxon-test for real data!")
            }
        }
        
    } else if (length(dim(fmridata)) <= 3) {
        # 2d or 3d
        p_dim = utils::head(dim(fmridata), n = -1)
        time_span = utils::tail(dim(fmridata), n = 1)
        fmri_mat = matrix(fmridata, nrow = prod(p_dim), ncol = time_span)
        p_vec = c()
        for (i in 1:nrow(fmri_mat)) {
            if (typeof(fmri_mat[i, ]) != "complex" & method %in% c("t-test", "wilcoxon-test")) {
                p_val = fmri_p_val(fmridata = fmri_mat[i, ], stimulus_idx = stimulus_idx, test_type = method)
                if (is.na(p_val) == TRUE) {
                  p_val = 1
                }
            } else if (typeof(fmri_mat[i, ]) == "complex" & method %in% c("HotellingsT2", "Wilks-Lambda", "gLRT")) {
                p_val = tryCatch({
                  as.numeric(fmri_complex_p_val(fmridata = fmri_mat[i, ], stimulus_idx = stimulus_idx, method = method, ons = ons, dur = dur)$p.value)
                }, error = function(e) {
                  return(1)
                })  # end of try-catch
            } else {
                # wrong test method input
                if (typeof(fmridata) == "complex") {
                  stop("Invalid test type! You can only use HotellingsT2, Wilks-Lambda or gLRT for complex data!")
                } else {
                  stop("Invalid test type! You can only use t-test or wilcoxon-test for real data!")
                }
            }
            p_vec = c(p_vec, p_val)
        }
        dim(p_vec) = p_dim
        return(p_vec)
    }
    
    
    
    # if data is 4d Get the dimension of the fmridata
    dim1 = dim(fmridata)[1]
    dim2 = dim(fmridata)[2]
    dim3 = dim(fmridata)[3]
    time_span = dim(fmridata)[4]
    
    # Build a 3d array to store the p-values
    p_val_3d = array(1, dim = c(dim1, dim2, dim3))
    
    # Store the indices for activated records
    on_idx = stimulus_idx
    # Unless specified, off_idx would be the complement
    off_idx = setdiff(1:time_span, on_idx)
    if (is.null(rest_idex) == FALSE) {
        off_idx = rest_idex
    }
    # if non-complex data
    if (typeof(fmridata) != "complex" & method %in% c("t-test", "wilcoxon-test", "on_off_diff", "HRF")) {
        
        if (method %in% c("t-test", "wilcoxon-test")) {
            for (x in 1:dim1) {
                for (y in 1:dim2) {
                  for (z in 1:dim3) {
                    if (is.null(mask) == FALSE) {
                      if (mask[x, y, z] == 1) {
                        p_val_3d[x, y, z] = fmri_p_val(fmridata[x, y, z, ], stimulus_idx = stimulus_idx, test_type = method)
                      }
                    } else {
                      p_val_3d[x, y, z] = fmri_p_val(fmridata[x, y, z, ], stimulus_idx = stimulus_idx, test_type = method)
                    }
                  }
                }
                p_val_3d[x, , ][is.na(p_val_3d[x, , ])] = 1
            }
        } else if (method == "on_off_diff" & is.null(mask) == FALSE) {
            
            # Store volume data for each voxel inside the mask
            volume_df = data.frame(x = NA, y = NA, z = NA, volume = NA)
            for (i in 1:dim1) {
                for (j in 1:dim2) {
                  for (k in 1:dim3) {
                    if (mask[i, j, k] == 1) {
                      volume_df = rbind(volume_df, c(i, j, k, fmri_on_off_volume(fmridata, x = i, y = j, z = k)))
                    }
                  }
                }
            }
            # Calculate the p-value using empirical CDF
            volume_df = volume_df %>% na.omit() %>% mutate(x = as.integer(x), y = as.integer(y), z = as.integer(z)) %>% mutate(p_val = 1 - 
                ecdf(volume_df$volume)(volume))
            
            
            for (i in 1:nrow(volume_df)) {
                p_val_3d[as.integer(volume_df[i, 1]), as.integer(volume_df[i, 2]), as.integer(volume_df[i, 3])] = volume_df$p_val[i]
            }
            
        } else if (method == "HRF" & is.null(mask) == FALSE) {
            p_val_3d = fmri_hrf_p_val(fmridata, mask = mask, ons = ons, dur = dur)
            
        } else {
            stop("Invalid test type! You can only use t-test or wilcoxon-test or on_off_diff or HRF for real data!")
        }
        
        
    } else {
        # if data is complex valued
        if (!method %in% c("HotellingsT2", "Wilks-Lambda", "gLRT")) {
            stop("Invalid test type! You can only use HotellingsT2 or Wilks-Lambda or gLRT for complex data!")
        }
        
        for (i in 1:dim1) {
            for (j in 1:dim2) {
                for (k in 1:dim3) {
                  tryCatch({
                    p = 1
                      # If use inputs a mask
                    if (is.null(mask) == FALSE) {
                      if (mask[i, j, k] == 1) {
                        p = as.numeric(fmri_complex_p_val(fmridata[i, j, k, ], stimulus_idx = stimulus_idx, method = method, ons = ons, dur = dur)$p.value)
                      }
                      # If no mask is used
                    } else {
                        p = as.numeric(fmri_complex_p_val(fmridata[i, j, k, ], stimulus_idx = stimulus_idx, method = method, ons = ons, dur = dur)$p.value)
                    }
                    
                    # p = as.numeric( (ICSNP::HotellingsT2(Y12 ~ labels))$p.value )
                    
                  }, error = function(e) {
                  })  # end of try-catch
                  
                  p_val_3d[i, j, k] = p
                }
            }
            p_val_3d[i, , ][is.na(p_val_3d[i, , ])] = 1
        }
    }
    
    # Post-hoc processing
    if (is.null(fdr_corr) == FALSE) {
        dim(p_val_3d) = c(dim1 * dim2 * dim3)
        p_val_3d = p.adjust(p_val_3d, method = fdr_corr, n = length(p_val_3d))
        dim(p_val_3d) = c(dim1, dim2, dim3)
    }
    
    if (is.null(spatial_cluster.size) == FALSE & is.null(spatial_cluster.thr) == FALSE) {
        spatial_cluster.filter = cluster.threshold(1 - p_val_3d, level.thr = spatial_cluster.thr, size.thr = spatial_cluster.size)
        p_val_3d = 1 - ((1 - p_val_3d) * spatial_cluster.filter)
    }
    return(p_val_3d)
}
