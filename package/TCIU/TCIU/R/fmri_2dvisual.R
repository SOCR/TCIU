#' @title visualization of the 2D brain (axial, sagittal and coronal) with the activated areas
#' @description a visualization method, using \code{ggplot2} to draw the brain 
#' from axial, sagittal and coronal view with activated area identified by p-values
#'
#' @param pval a 3D array of p-values used to plot activated area of the brain
#' @param axis_ls a list with two elements. The first element is the character of 'x', 'y', 'z'.
#' The second element is an integer showing a specific slice on the fixed axis identified in the first element.
#' @param hemody_data a parameter to have the plot with/without hemodynamic contour. The default is NULL to make the plot 
#' without hemodynamic contour, otherwise assign a 3D array of the hemodynamic data.
#' @param mask a 3D nifti or 3D array of data to show the shell of the brain
#' @param p_threshold NULL or a numeric value that can be selected randomly below 0.05 to 
#' drop all p-values above the threshold. If 'low5_percent' method is used, 
#' make 'p_threshold' as NULL. The default is 0.05.
#' @param legend_show a logical parameter to specify whether the final plot has legend
#' @param method a string that represents method for the plot. 
#' There are 3 options: 'min_max', 'scale_p' and 'low5_percent'. The default is 'scale_p'.
#' 'min_max' is to draw plot based on the color scale of the minimum and maximum of the p value; 
#' 'scale_p' is to draw the plot with fixed color scale for fixed range of p value; 
#' 'low5_percent' is to draw the plot for the smallest 5 percent of p value when all the p values are not significant.
#' @param color_pal the name of the color palettes provided by \code{RColorBrewer}. The default is "YlOrRd".
#' @param multi_pranges an option under 'scale_p' method to decide whether there are at most 9 colors 
#' in the legend for the ranges of p value, or at most 4 colors. 
#' The default is TRUE, choosing the larger number of colors for the plot.
#' @param mask_width a numeric value to specify the width of mask contour. The default is 1.5.
#' 
#' @details
#' The function \code{fmri_2dvisual} is used to find activated part of the brain 
#' based on given p values from sagittal, axial and coronal view. When providing input of 
#' the p-values, the specific plane and index to slice on, the mask data and 
#' the hemodynamic data of the brain, a plot will be generated with the heat map 
#' for the activated parts, the black contour showing the position of the brain, 
#' and the blue contour representing the hemodynamic contour.
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return a plot drawn by \code{ggplot2}
#'
#' @examples
#' # sample 3D data of mask provided by the package
#' dim(mask)
#' # sample 3D p value provided by the package
#' dim(phase2_pval)
#' 
#' # plot the sagittal, coronal and axial view of this p value generated from the brain fMRI data
#' \donttest{
#' fmri_2dvisual(phase2_pval, list('x',40), hemody_data=NULL, mask=mask, p_threshold=0.05)
#' }
#' 
#' @export
#'
#' @import RColorBrewer fancycut scales ggpubr dplyr
#' @importFrom pracma ceil
#' @importFrom ggplot2 coord_fixed theme geom_contour aes ggtitle scale_fill_gradient scale_alpha geom_tile ggplot labs geom_col element_text


fmri_2dvisual = function(pval,
                         axis_ls,
                         hemody_data = NULL, 
                         mask,
                         p_threshold = 0.05,
                         legend_show = TRUE,
                         method = "scale_p",
                         color_pal = "YlOrRd",
                         multi_pranges = TRUE,
                         mask_width = 1.5){
    
    floor_dec = function(x, level=1) round(x - 5*10^(-level-1), level)
    ceiling_dec = function(x, level=1) round(x + 5*10^(-level-1), level)
    # "prange_extract()" is a small function, specialized to generate the unique name of p values intervals
    # based on the input of the unique name of 1-p values intervals
    # the input of this function are the characters of 1 - p values range intervals, 
    # the output are the characters of p values range intervals
    prange_extract = function(name_range_one_minus_p){
      
      numeric_range_one_minus_p = regmatches(name_range_one_minus_p, 
                                             gregexpr("[[:digit:]\\.]+", 
                                                      name_range_one_minus_p))
      p_range = c()
      for(i in 1:length(numeric_range_one_minus_p)){
        p_2range = 1-as.numeric(numeric_range_one_minus_p[[i]])
        p_range[i] = sprintf("(%.1e, %.1e]", p_2range[2], p_2range[1])
      }
      return(p_range)
    }
    
    
    # generate customized error message based on different wrong input
    if((class(pval)!="array") |
       (length(dim(pval))!=3)){
      stop("'pval' should be a 3d array.")
    }else if((class(axis_ls)!="list" | 
              ((axis_ls[[1]] %in% c("x", "y", "z")) == FALSE))){
      stop("'axis_ls' should be a list, with the first element as a string from 'x', 'y' or 'z'.")
    }else if(is.null(hemody_data) != TRUE){
      if((class(hemody_data)!="array") | 
         (length(dim(hemody_data))!=3)){
        stop("If 'hemody_data' is not NULL, then it should be a 3d array.")
      }
    }else if(((class(mask) %in% c("array", "nifti")) != TRUE) |
             (length(dim(mask))!=3)){
      stop("'mask' should be a 3d array.")
    }else if(is.null(p_threshold) != TRUE){
      if((class(p_threshold) != "numeric") | 
         (p_threshold > 0.05) | (p_threshold <= 0)){
        stop("'p_threshold should be NULL or a numeric value in range of (0, 0.05].'")
      }
    }else if(class(legend_show) != "logical"){
      stop("'legend_show' should be a logical TRUE or FALSE.")
    }else if((method %in% c("scale_p", "min_max", "low5_percent")) != TRUE){
      stop("'method' should only choose from 'scale_p', 'min_max' or 'low5_percent'.")
    }
    
    one_minus_p_3d = 1 - pval
    
    dim_xyz = dim(one_minus_p_3d)
    dim_x = dim_xyz[1]
    dim_y = dim_xyz[2]
    dim_z = dim_xyz[3] 
    
    dim_i_bound = switch(axis_ls[[1]], 
                         "x" = {dim_x},
                         "y" = {dim_y},
                         "z" = {dim_z})
    if((ceil(axis_ls[[2]]) != axis_ls[[2]]) | 
       (axis_ls[[2]]>dim_i_bound) | (axis_ls[[2]]<=0)){
      stop("'The second element of axis_ls should be an integer that is not out of range.'")
    }

    
    # based on given x|y|z, compute and manage needed data for the plot correspondingly
    switch(axis_ls[[1]],
           "x" = {
             # rehsape 3 dimension 1-p data to have it saved in dataframe
             dim(one_minus_p_3d) = c(dim_x, dim_y*dim_z)
             # this part manage 1-p data for motor area heatmap
             one_minus_p_df = mask_df =
               expand.grid(x = seq(1, dim_y, by = 1), 
                           y = seq(1, dim_z, by = 1))
             one_minus_p_df$one_minus_p = one_minus_p_3d[axis_ls[[2]],]
             # this part manage mask data for brain shell contour
             contour_mask_2d = mask[axis_ls[[2]],,]
             # this part manage mask data for hemodynamic contour
             contour_hemody_2d = hemody_data[axis_ls[[2]],,]
             contour_hemody_df = expand.grid(x=seq(1, dim_y, by = 1), 
                                             y=seq(1, dim_z, by = 1))
             # title and axis name for future plot
             title = paste0('Sagittal View of Brain for x=', axis_ls[[2]])
             plot_lab = labs(x="y", y="z", fill="p value")
           },
           
           "y" = {
             # rehsape 3 dimension 1-p data to have it saved in dataframe
             one_minus_p_3d = aperm(one_minus_p_3d, c(1,3,2))
             dim(one_minus_p_3d) = c(dim_x*dim_z, dim_y)
             # this part manage 1-p data for motor area heatmap
             one_minus_p_df = mask_df =
               expand.grid(x = seq(1, dim_x, by = 1), 
                           y = seq(1, dim_z, by = 1))
             one_minus_p_df$one_minus_p = one_minus_p_3d[, axis_ls[[2]]]
             # this part manage mask data for brain shell contour
             contour_mask_2d = mask[, axis_ls[[2]],]
             # this part manage mask data for hemodynamic contour
             contour_hemody_2d = hemody_data[,axis_ls[[2]],]
             contour_hemody_df = expand.grid(x=seq(1, dim_x, by = 1), 
                                             y=seq(1, dim_z, by = 1))
             
             title = paste0("Coronal View of Brain for y=", axis_ls[[2]])
             plot_lab = labs(x="x", y="z", fill="p value")
           },
           
           "z" = {
             # rehsape 3 dimension 1-p data to have it saved in dataframe
             dim(one_minus_p_3d) = c(dim_x*dim_y, dim_z)
             # this part manage 1-p data for motor area heatmap
             one_minus_p_df = mask_df =
               expand.grid(x = seq(1, dim_x, by = 1), 
                           y = seq(1, dim_y, by = 1) )
             one_minus_p_df$one_minus_p = one_minus_p_3d[, axis_ls[[2]]]
             # this part manage mask data for brain shell contour
             contour_mask_2d = mask[,, axis_ls[[2]]]
             # this part manage mask data for hemodynamic contour
             contour_hemody_2d = hemody_data[,,axis_ls[[2]]]
             contour_hemody_df = expand.grid(x=seq(1, dim_x, by = 1), 
                                             y=seq(1, dim_y, by = 1))
             # title and axis name for future plot
             title = paste0('Axial View of Brain for z=', axis_ls[[2]])
             plot_lab = labs(x="x", y="y", fill="p value")
             
           }
    )
    
    
    # put mask data into 1-p dataframe according to the x, y value
    # so that at specific voxel the values for 1-p and mask value are fixed
    mask_df$maskval = NA
    idx = 1
    for(j in unique(mask_df$y)){
      for(i in unique(mask_df$x)){
        mask_df$maskval[idx] = contour_mask_2d[i, j]
        idx = idx + 1
      }
    }
    
    contour_bin = 1

    # based on given method, make the integrated plot correspondingly
    if (method %in% c("scale_p", "low5_percent")){
      pval_df = fmri_3dvisual(pval, mask, p_threshold, 
                            method, color_pal = color_pal, 
                            multi_pranges)$pval_df
  
      switch(method,
             "scale_p" = {
               # drop the corresponding rows if p value is larger than the given p_threshold
               one_minus_p_df =  dplyr::filter(one_minus_p_df, one_minus_p_df$one_minus_p >= 1 - p_threshold)
               one_minus_p_df$p_val = 1 - one_minus_p_df$one_minus_p
               if(nrow(one_minus_p_df) == 0){
                 stop("There does not exist any valid data under this 'p_threshold'. Please change a 'p_threshold' or change the method to 'low5_percent' and make 'p_threshold as NULL'.")
               }
               one_minus_p_df$cut_invs = wafflecut(one_minus_p_df$p_val, 
                                                   c('[0,1e-8]', '(1e-8, 1e-7]', '(1e-7, 1e-6]', '(1e-6, 1e-5]', 
                                                     '(1e-5, 1e-4]', '(1e-4, 1e-3]', '(1e-3, 1e-2]', '(1e-2, 5e-2]'))
               one_minus_p_df$colorgrp = as.numeric(one_minus_p_df$cut_invs)
               if(multi_pranges == FALSE){
                 one_minus_p_df$cut_invs = wafflecut(one_minus_p_df$p_val, 
                                                     c('[0,1e-7]', '(1e-7, 1e-5]', 
                                                       '(1e-5, 1e-3]', '(1e-3, 5e-2]'))
                 one_minus_p_df$colorgrp= as.numeric(one_minus_p_df$cut_invs)
               }
             },
             
             "low5_percent" = {
               
               one_minus_p_df$p_val = 1 - one_minus_p_df$one_minus_p
               quantile5_pt = quantile(one_minus_p_df$p_val, 0.05)
               one_minus_p_df = dplyr::filter(one_minus_p_df, one_minus_p_df$p_val <= quantile5_pt)
               p_val_max_minus_min = max(one_minus_p_df$p_val)-min(one_minus_p_df$p_val)
               if (p_val_max_minus_min != 0){
                 cut_pts = min(one_minus_p_df$p_val) + 0:9 * (p_val_max_minus_min) / 9
                 cut_pts[length(cut_pts)] = ceiling_dec(cut_pts[length(cut_pts)], 2)
                 cut_pts[1] = floor_dec(cut_pts[1], 2)
                 # add lapply to cut pts[2:-1] to round 3
                 cut_pts = unique(sapply(cut_pts, function(x) round(x, 3)))
                 invs_total = levels(cut(one_minus_p_df$p_val, cut_pts, right=TRUE))
                 invs_total[1] = sub("\\(", "[", invs_total[1])
                 one_minus_p_df$cut_invs = wafflecut(one_minus_p_df$p_val, invs_total)
                 one_minus_p_df$colorgrp = as.numeric(one_minus_p_df$cut_invs)##
               }else{
                 one_minus_p_df$cut_invs = as.character(round(min(one_minus_p_df$p_val), 3))
                 for (i in unique(pval_df$cut_invs)){
                   if(length(levels(wafflecut(as.numeric(unique(one_minus_p_df$cut_invs)),
                                              c(as.character(i))))) == 1){
                     one_minus_p_df$colorgrp = unique(pval_df[which(pval_df$cut_invs==i),]$colorgrp)
                   }
                 }
                 
               }
             }
             
      )
      one_minus_p_df = merge(x = one_minus_p_df,
                             y = unique(pval_df%>%dplyr::select(colorgrp, corresp_color)),
                             by = "colorgrp", all.x = TRUE)
      
      motorinbrain = 
        ggplot() +
        geom_tile(one_minus_p_df, mapping=aes(x=x, y=y, fill=corresp_color)) + 
        scale_fill_identity("p value", 
                            labels= rev(as.character(unique(one_minus_p_df$cut_invs))), 
                            breaks = rev(unique(one_minus_p_df%>%dplyr::select(corresp_color))[[1]]), 
                            guide = "legend") +
        plot_lab + 
        geom_contour(mask_df, mapping = aes(x=mask_df$x, y=mask_df$y, z=mask_df$maskval), 
                     colour = "black", bins = contour_bin, size=mask_width) + 
        ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
      
    }else if (method == "min_max"){
      
      # drop the corresponding rows if p value is larger than the given p_threshold
      one_minus_p_df =  dplyr::filter(one_minus_p_df, one_minus_p_df$one_minus_p >= 1 - p_threshold)
      one_minus_p_df$p_val = 1 - one_minus_p_df$one_minus_p
      if(nrow(one_minus_p_df) == 0){
        stop("There does not exist any valid data under this 'p_threshold'. Please change a 'p_threshold' or change the method to 'low5_percent' and make 'p_threshold as NULL'.")
      }
      
      idx = one_minus_p_df$one_minus_p != 0
      one_minus_p_min = min(one_minus_p_df$one_minus_p[idx])
      one_minus_p_max = max(one_minus_p_df$one_minus_p[idx])
      break_pt = c(one_minus_p_min, (one_minus_p_min+one_minus_p_max)/2, one_minus_p_max)
      motorinbrain = 
        ggplot() + 
        geom_tile(one_minus_p_df, mapping = aes(x=x, y=y, fill=one_minus_p, alpha=one_minus_p)) + 
        scale_alpha(range = c(0, 1), guide=F) +
        plot_lab +
        scale_fill_gradient(low = "yellow", high = "red", 
                            limits=c(break_pt[1], break_pt[length(break_pt)]), 
                            oob=squish, breaks=break_pt, 
                            labels=format(round(1-break_pt, 4), nsmall = 4)) +
        geom_contour(mask_df, mapping = aes(x=mask_df$x, y=mask_df$y, z=mask_df$maskval), 
                     colour = "black", bins = contour_bin, size=mask_width) + 
        ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    }
    
    # manage hemody data to make it dataframe for the plot
    if(is.null(hemody_data) == FALSE){
      # unsmooth data for blood flow in brain
      contour_hemody_df$modval = NA
      idx = 1
      for(j in unique(contour_hemody_df$y)){
        for(i in unique(contour_hemody_df$x)){
          contour_hemody_df$modval[idx] = contour_hemody_2d[i, j]
          idx = idx + 1
        }
      }
      
      motorinbrain = motorinbrain +
        geom_contour(contour_hemody_df, mapping = aes(x=contour_hemody_df$x, y=contour_hemody_df$y, z=contour_hemody_df$modval), 
                     colour = "blue4", bins = 8, size=0.7)
    }
    
    # decide whether the final plot include the legend or not
    if(legend_show){
      return(motorinbrain + coord_fixed())
             get_legend(motorinbrain)
    }else{
      return(motorinbrain + coord_fixed() + 
             theme(legend.position="none"))
    }
    
}
