#' @title comparison between 3d visualization for p-values
#' @description a visualization method, use \code{plotly} to compare 
#' the activated parts inside the brain, using two sets of color palettes. 
#' The activated parts are localized with different p values.
#'
#' @param pval_3d_ls a list of two element, each element is a 3D array of p-values used to plot activated area of the brain
#' @param mask a 3D nifti or 3D array of data to show the shell of the brain
#' @param p_threshold NULL or a numeric value that can be selected randomly below 0.05 to 
#' drop insignificant p-values of no need or drop no p-values. If 'low5_percent' method is used, 
#' make 'p_threshold' as NULL. The default is 0.05.
#' @param method_ls a string that represents method for the plot. 
#' There are 2 options: 'scale_p' and 'low5_percent'. The default is 'scale_p'. 
#' 'scale_p' is to draw the plot with fixed color scale for fixed range of p value; 
#' 'low5_percent' is to draw the plot for the smallest 5 percent of p value 
#' when all the p values are not significant
#' @param color_pal_ls a list of two element. Each element is the name of 
#' the color palettes provided by \code{RColorBrewer}. The default is list('YlOrRd', 'YlGnBu').
#' @param multi_pranges an option under 'scale_p' method to decide whether there are at most 9 colors 
#' in the legend for the ranges of p value, or at most 4 colors. 
#' The default is TRUE, choosing the larger number of colors for the plot.
#' 
#' @details
#' The function \code{fmri_pval_comparison_3d} is used to visualize and compare the 3D plots of 
#' the activated parts in one brain shell. The activated parts are plotted 
#' based on p-values provided. Note that this comparison can only be made when the masks 
#' of the two p values are the same. When providing input of two set of the 3D array of p-values, 
#' corresponding p threshold for each p value data, and the method to draw the plot, 
#' the plot will be generated with one brain shell and two groups of activated parts 
#' in two sets of color palettes. The size and color of the scatter points represent 
#' different stimulated levels of the activated parts.
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return a plot drawn by \code{plotly}
#'
#' @examples
#' # sample 3D data of mask provided by the package
#' dim(mask)
#' # sample 3D p value provided by the package
#' dim(phase2_pval)
#' dim(phase3_pval)
#' 
#' fmri_pval_comparison_3d(list(phase2_pval, phase3_pval), mask, 
#'                         list(0.05, 0.05), list('scale_p', 'scale_p'), multi_pranges=FALSE)
#' @export
#'
#' @import plotly dplyr scales RColorBrewer fancycut geometry

fmri_pval_comparison_3d = function(pval_3d_ls,
                                   mask,
                                   p_threshold,
                                   method_ls,
                                   color_pal_ls = list("YlOrRd", "YlGnBu"),
                                   multi_pranges = TRUE) {
    
    plot_comp1 = fmri_3dvisual(pval_3d_ls[[1]], mask, p_threshold[[1]], method_ls[[1]], color_pal = color_pal_ls[[1]], multi_pranges)$plot
    df_pval2 = fmri_3dvisual(pval_3d_ls[[2]], mask, p_threshold[[2]], method_ls[[2]], color_pal = color_pal_ls[[2]], multi_pranges)$pval_df
    
    if (multi_pranges == FALSE) {
        color_choice = rev(brewer.pal(9, color_pal_ls[[2]]))[c(FALSE, TRUE)]
    } else {
        color_choice = rev(brewer.pal(9, color_pal_ls[[2]]))
    }
    
    for (i in sort(unique(df_pval2$colorgrp))) {
        
        pts_grp = df_pval2[df_pval2$colorgrp == i, ]
        plot_comp1 = add_markers(plot_comp1, data = pts_grp, x = ~x, y = ~y, z = ~z, mode = "markers", name = paste0("p value in ", 
            df_pval2$cut_invs)[i], marker = list(opacity = 0.6, symbol = "triangle", size = rev(c(1:length(color_choice)))[i] + 
            2 * as.numeric(I(multi_pranges == FALSE)), color = ~corresp_color, line = list(color = ~corresp_color, width = 0.3)))
    }
    
    return(plot_comp1)
    
}
