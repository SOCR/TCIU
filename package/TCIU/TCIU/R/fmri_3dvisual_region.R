#' @title visualization of the 3D brain with the activated areas by regions
#' @description an improved visualization method of \code{fmri_3dvisual}, using \code{plotly} to 
#' draw the 3D plot of the brain with the activated areas region by region
#'
#' @param pval a 3D or 1D or a list of two 3D array of p-values used to plot activated area of the brain
#' @param mask a 3D nifti or 3D array of data to show the regions of the brain
#' @param label_index a 1D array listing  the label number in the mask
#' @param label_name a 1D array corresponding to the name of the label number in the mask
#' @param top_num NULL or a numeric value that used for 1D p-values. If specified, the output will 
#' show the top num significant regions. The default is NULL.
#' @param p_threshold NULL or a numeric value that used for 3D p-values can be selected randomly below 0.05 to 
#' drop all p-values above the threshold. If 'low5_percent' method is used, 
#' make 'p_threshold' as NULL. The default is 0.05.
#' @param method a string that represents method for the 3D p-values plot.
#' There are 2 options: 'scale_p' and 'low5_percent'. The default is 'scale_p'. 
#' 'scale_p' is to draw the plot with fixed color scale for fixed range of p value.
#' 'low5_percent' is to draw the plot for the smallest 5 percent of p value when all the p values are not significant.
#' @param multi_pranges an option under 'scale_p' method to decide whether there are at most 9 colors 
#' in the legend for the ranges of 3D p-values, or at most 4 colors. 
#' The default is TRUE, choosing the larger number of colors for the plot.
#' @param color_pal the name of the color palettes provided by \code{RColorBrewer}. The default is "YlOrRd".
#' @param rank the method that how the trace is ranked. The default is NULL.
#' There are 2 options: 'value' and a vector.
#' 'value' is to draw the 1D p-values by the values from smallest to largest.
#' a vector is to specific the rank of the regions in 3D p-values plot.
#' @param title the title of the plot. The default is NULL.
#' 
#' @details
#' The function \code{fmri_3dvisual_region} is used to visualize the 3D plot of the brain 
#' with activated parts region by region. When providing a 1D/3D p-values data, a 3D interactive
#' plot with surface of the brain shell will be generated with either scatter points representing 
#' different stimulated levels or large color pieces representing different regions of the brain. 
#' When providing a list of two 3D array of p-values, two 3D interactive brains with different scatter
#' points corresponding to the two input 3D p-values will be given.
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return the 3d plot of the fMRI data drawn by \code{plotly}
#' 
#' @examples
#' # sample label vector provided in the package
#' label_index = mask_dict$index
#' label_name = as.character(mask_dict$name)
#' label_mask = mask_label
#' \donttest{
#' fmri_3dvisual_region(phase1_pval, label_mask, label_index,
#'                      label_name, title = "phase1 p-values")
#' fmri_3dvisual_region(phase1_pval, label_mask, label_index,
#'                      label_name, 5, title = "phase1 top five p-values", rank = "value")
#' 
#' # for 3D visualization, user needs to include empty region in the label
#' label_index = c(0, label_index)
#' label_name = c("empty", label_name)
#' fmri_3dvisual_region(phase2_pval, label_mask, label_index,
#'                      label_name, title = "phase2 p-values", rank = c(1:length(label_name)))
#' fmri_3dvisual_region(list(phase2_pval,phase3_pval), label_mask, label_index,
#'                      label_name, title = "phase2&3 p-values")
#' }
#' @export
#'
#' @import plotly scales RColorBrewer fancycut geometry dplyr
#' @importFrom grDevices colorRampPalette contourLines
#' 
#' 

fmri_3dvisual_region = function(pval,
                                mask,
                                label_index,
                                label_name,
                                top_num = NULL,
                                p_threshold = 0.05,
                                method = "scale_p",
                                multi_pranges = TRUE,
                                color_pal = "YlOrRd",
                                rank = NULL,
                                title = NULL){
  floor_dec = function(x, level=1) round(x - 5*10^(-level-1), level)
  ceiling_dec = function(x, level=1) round(x + 5*10^(-level-1), level)
  # get the dim of the mask
  xdim = dim(mask)[1]
  ydim = dim(mask)[2]
  zdim = dim(mask)[3]  
  fig = plot_ly()

  if (!is.null(title)){fig = fig%>%layout(title = title)}
  flag = FALSE
  if (is.list(pval)){
    pval1 = pval[[1]]
    pval2 = pval[[2]]
    flag = TRUE
    pval = array(1, dim = c(2*xdim,ydim,zdim))
    for (i in 1:xdim){
      for (j in 1:ydim){
        for (k in 1:zdim){
            pval[i,j,k] = pval1[i,j,k]
            pval[i+xdim,j,k] = pval2[i,j,k]
          }
        }
      }
  }
  # "newContour()" is a small function, specialized to generate boundary point for the 2d z-plane 
  # of a 3d structure based on the the location of z
  # the input of this function are the location of z and the 3d structure as a 3d array 
  # the output is a dataframe which contains the index of the boundary points for 
  # 2d plane of the 3d structure at location z
  newContour = function(z, struc_3d = newMask){
    dim_maskz = dim(struc_3d[,,z])
    x = seq(1, dim_maskz[1])
    y = seq(1, dim_maskz[2])
    contour_pt = contourLines(x = x, y = y, z = struc_3d[,,z])
    if (length(contour_pt) != 0){
      boundry_pt_df = cbind(contour_pt[[1]]$x, contour_pt[[1]]$y, z)
      return(boundry_pt_df)
    }
  }
  
  newMask = array(0, dim = c(xdim,ydim,zdim))
  # add the outer mask shell to the output
  for (i in 1:xdim){
    for (j in 1:ydim){
      for (k in 1:zdim){
        if (mask[i,j,k] != 0){
          newMask[i,j,k] = 1
        }
      }
    }
  }
  outerBoundry = data.frame(do.call(rbind, lapply(1:dim(newMask)[3], newContour)))
  colnames(outerBoundry) = c("x", "y", "z")
  outerTri = delaunayn(outerBoundry)
  facecolor = rep("#F0FFFF",length(outerTri[, 1]))
  fig = fig %>% 
    add_trace(outerBoundry, 
              x = outerBoundry[, 1], y = outerBoundry[, 2], z = outerBoundry[, 3],
              i = outerTri[, 1]-1, j = outerTri[, 2]-1, k = outerTri[, 3]-1,
              type = "mesh3d", opacity = 0.01, name = "Brain Shell",
              #facecolor = facecolor,
              contour = list(show = TRUE, color="#000", width=15))
              #showlegend = TRUE)
  if (flag){
    newMask = array(0, dim = c(2*xdim,ydim,zdim))
    # add the outer mask shell to the output
    for (i in 1:xdim){
      for (j in 1:ydim){
        for (k in 1:zdim){
          if (mask[i,j,k] != 0){
            newMask[i+xdim,j,k] = 1
          }
        }
      }
    }
    outerBoundry = data.frame(do.call(rbind, lapply(1:dim(newMask)[3], newContour)))
    colnames(outerBoundry) = c("x", "y", "z")
    outerTri = delaunayn(outerBoundry)
    facecolor = rep("#F0FFFF",length(outerTri[, 1]))
    fig = fig %>% 
      add_trace(outerBoundry, 
                x = outerBoundry[, 1], y = outerBoundry[, 2], z = outerBoundry[, 3],
                i = outerTri[, 1]-1, j = outerTri[, 2]-1, k = outerTri[, 3]-1,
                type = "mesh3d", opacity = 0.01, name = "Brain Shell",
                #facecolor = facecolor,
                contour = list(show = TRUE, color="#000", width=15))
                #showlegend = TRUE)
  }
  if (flag){
    innerMask = array(0, dim = c(2*xdim,ydim,zdim))
    # change the inner mask number to begin with 1 
    for (i in 1:xdim){
      for (j in 1:ydim){
        for (k in 1:zdim){
          temp = mask[i,j,k]
          for (index in 1:length(label_index)){
            if (temp == label_index[index]){
              innerMask[i,j,k] = index
              innerMask[i+xdim,j,k] = index
              break
            }
          }
        }
      }
    }
  } else {
    innerMask = array(0, dim = c(xdim,ydim,zdim))
    # change the inner mask number to begin with 1 
    for (i in 1:xdim){
      for (j in 1:ydim){
        for (k in 1:zdim){
          temp = mask[i,j,k]
          for (index in 1:length(label_index)){
            if (temp == label_index[index]){
              innerMask[i,j,k] = index
              break
            }
          }
        }
      }
    }
  }
  if (is.null(dim(pval))){
    # 1D p-value
    # generate length(pval) different colors
    # color_choice[1] is the smallest pval
    color_choice = colorRampPalette(brewer.pal(9,color_pal))(length(pval))
    if (!is.null(top_num)){
      color_choice = colorRampPalette(brewer.pal(9,color_pal))(top_num)
    }
    color_choice = rev(color_choice)
    # sort the color index
    oldPval = pval
    newPval = sort(oldPval)
    colorPval = array(0,dim=length(pval))

    ranking = array(0,dim = length(pval))
    for (i in 1:length(pval)){
      for (j in 1:length(pval)){
        if (oldPval[i] == newPval[j]){
          if (!is.null(top_num) && j>top_num){
            colorPval[i] = -1
          } else {
            colorPval[i] = j
          }
          ranking[j] = i
        }
      }
    }
    
    if (!is.null(rank) && rank == "value"){
      ranking = ranking
    } else {
      ranking = c(1:length(pval))
    }

    for (x in ranking){
      if (colorPval[x] > 0){
        newMask = array(0,dim=c(xdim, ydim, zdim))
        # now we get new mask with only "L superior frontal gyrus" to be 1
        for (i in 1:xdim){
          for (j in 1:ydim){
            for (k in 1:zdim){
              if (innerMask[i,j,k]==x){
                newMask[i,j,k] = 1
              }
            }
          }
        }
        # then we can calculate the contour
        boundryContour = data.frame(do.call(rbind, lapply(1:dim(newMask)[3], newContour)))
        colnames(boundryContour) = c("x", "y", "z")
        triContour = delaunayn(boundryContour)
        faceColor = rep(color_choice[colorPval[x]],length(triContour[, 1]))
        fig = fig %>%
          add_trace(boundryContour, 
                    x=boundryContour[, 1], y=boundryContour[, 2], z=boundryContour[, 3],
                    i = triContour[, 1]-1, j = triContour[, 2]-1, k = triContour[, 3]-1,
                    type = "mesh3d", 
                    opacity = 1,
                    name = label_name[x],
                    hovertemplate  = paste("p-value:",pval[x]),
                    facecolor = faceColor,
                    contour = list(show = TRUE, color="#000", width=15),
                    showlegend=TRUE) 
      }
    }
    
  } else {
    # 3D p-value
    # this part is to process p-value data to add the points of the motor area in the mesh brain plot
    # this plot can only be made on 1 - p value
    dim_pval = dim(pval)
    # convert 3d array 1 - p value to a dataframe to be the input for plotly
    pval_df = 
      data.frame(x = rep(c(1:dim_pval[1]), dim_pval[2]*dim_pval[3]),
                 y = rep(rep(c(1:dim_pval[2]), each=dim_pval[1]), dim_pval[3]),
                 z = rep(c(1:dim_pval[3]), each=dim_pval[1]*dim_pval[2]))
    names(pval_df) = c("x", "y", "z")
    pval_df$p_val = as.vector(pval)
    
    switch (method,
            "scale_p" = {
              pval_df = dplyr::filter(pval_df, pval_df$p_val <= p_threshold)
              
              if(multi_pranges == FALSE){
                pval_df$cut_invs = wafflecut(pval_df$p_val, 
                                             c('[0,1e-7]', '(1e-7, 1e-5]', '(1e-5, 1e-3]', '(1e-3, 5e-2]'))
              } else {
                pval_df$cut_invs = wafflecut(pval_df$p_val, 
                                             c('[0,1e-8]', '(1e-8, 1e-7]', '(1e-7, 1e-6]', '(1e-6, 1e-5]', 
                                               '(1e-5, 1e-4]', '(1e-4, 1e-3]', '(1e-3, 1e-2]', '(1e-2, 5e-2]'))
              }
              # consider whether we can decide to change the interval cut based on the condition of p-value provided
              pval_df$colorgrp = as.numeric(pval_df$cut_invs)
              if (length(unique(pval_df$colorgrp)) <= 4){
                print("The number of p ranges is originally less than 4, 
                       it is recommended to turn 'multi_pranges' as FALSE.")
              }
            },
            
            "low5_percent" = {
              quantile5_pt = quantile(pval_df$p_val, 0.05)
              pval_df = dplyr::filter(pval_df, pval_df$p_val <= quantile5_pt)
              p_val_max_minus_min = max(pval_df$p_val)-min(pval_df$p_val)
              if (p_val_max_minus_min != 0){
                cut_pts = min(pval_df$p_val) + 0:9 * (p_val_max_minus_min) / 9
                cut_pts[length(cut_pts)] = ceiling_dec(cut_pts[length(cut_pts)], 2)
                cut_pts[1] = floor_dec(cut_pts[1], 2)
                # add lapply to cut pts[2:-1] to round 3
                cut_pts = unique(sapply(cut_pts, function(x) round(x, 3)))
                invs_total = levels(cut(pval_df$p_val, cut_pts, right=TRUE))
                invs_total[1] = sub("\\(", "[", invs_total[1])
                pval_df$cut_invs = wafflecut(pval_df$p_val, invs_total)
                pval_df$colorgrp = as.numeric(pval_df$cut_invs)
              }else{
                pval_df$cut_invs = as.character(round(min(pval_df$p_val), 3))
                pval_df$colorgrp = 1
              }
              
            }
    )
    
    if (multi_pranges == FALSE){
      color_choice = rev(brewer.pal(9, color_pal))[c(FALSE, TRUE)]
    }else{
      color_choice = rev(brewer.pal(9, color_pal))
    }
    
    pval_df$corresp_color = color_choice[pval_df$colorgrp]
    pvalX = pval_df$x
    pvalY = pval_df$y
    pvalZ = pval_df$z
    groupname = c()
    groupindex = c()
    for (i in 1:length(pvalX)){
      tempIndex = innerMask[pvalX[i],pvalY[i],pvalZ[i]]
      if (!is.element(tempIndex, groupindex)){
        groupindex = c(groupindex, tempIndex)
      }
      add_name =  label_name[tempIndex]
      groupname = c(groupname, add_name)
    }
    pval_df$groupname = groupname
    for (name in sort(unique(pval_df$groupname))){
      pts_grp = pval_df[pval_df$groupname == name, ]
      # plot the small p-values first
      pts_grp = pts_grp[sort(pts_grp$colorgrp,index.return=TRUE)$ix,]
      size = c()
      for (i in 1:length(pts_grp$colorgrp)){
        size = c(size, rev(c(1:length(color_choice)))[pts_grp$colorgrp[i]]+2*as.numeric(I(multi_pranges==FALSE))+2)
      }
      
      fig = add_markers(fig, data=pts_grp,
                        x=~x, y=~y, z=~z, mode="markers", 
                        legendgroup = name,
                        marker = list(opacity=0.6, symbol="circle", 
                                      size=size,
                                      color=~corresp_color,
                                      line = list(color = ~corresp_color,
                                                  width = 0.3)),
                        showlegend = FALSE
      )
      
    }
    
    
    # add the corresponding inner mask shell
    if (!is.null(rank)){
      ranking = c()
      for (i in rank){
        if (is.element(i, groupindex)){
          ranking = c(ranking, i)
        }
      }
    } else {
      ranking = groupindex
    }
    for (x in ranking){
      newMask = array(0,dim=c(xdim, ydim, zdim))
      # now we get new mask with only "L superior frontal gyrus" to be 1
      for (i in 1:xdim){
        for (j in 1:ydim){
          for (k in 1:zdim){
            if (innerMask[i,j,k]==x){
              newMask[i,j,k] = 1
              break
            }
          }
        }
      }
      # then we can calculate the contour
      boundryContour = data.frame(do.call(rbind, lapply(1:dim(newMask)[3], newContour)))
      colnames(boundryContour) = c("x", "y", "z")
      triContour = delaunayn(boundryContour)
      faceColor = rep("#F0FFFF",length(triContour[, 1]))
      fig = fig %>%
        add_trace(boundryContour, 
                  x=boundryContour[, 1], y=boundryContour[, 2], z=boundryContour[, 3],
                  i = triContour[, 1]-1, j = triContour[, 2]-1, k = triContour[, 3]-1,
                  type = "mesh3d", 
                  opacity = 0,
                  facecolor = faceColor,
                  name = label_name[x],
                  legendgroup = label_name[x],
                  contour = list(show = TRUE, color="#000", width=15),
                  showlegend=TRUE) 
    }
  }
  
  return (fig)
}
