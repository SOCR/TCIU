#' @title post-hoc process for p values
#' @description This function is used to conduct the post-hoc process (i.e. FDR correction and spatial clustering) for a 3-dimensional p-value array.
#'
#' @param p_val_3d an array which contains the p-values as the result of fMRI statistical tests.
#' @param fdr_corr The default is NULL. Input 'fdr' to conduct FDR correction.
#' @param spatial_cluster.thr The default is NULL. Together with spatial_cluster.size are used to filter contiguous clusters of locations in a 3D array that are below some threshold and with some minimum size.
#' @param spatial_cluster.size The default is NULL. The size of spatial cluster.
#' @param show_comparison The default is FALSE. If TRUE, the output would display the comparison between raw and processed p-values.
#' @param ...  One can specify breaks etc. to modify the comparison histogram in ggplot2.

#' @details
#' The function \code{fmri_post_hoc} would help do the FDR correction and spatial clustering for a 3d p-value array. The FDR correction controls for a low proportion of false positives, while the spatial clustering part help filter out all sparse p-values that are not in specified clusters.
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return 3D p-values after FDR correction or spatial clustering
#' @export
#' 
#' @import tidyr
#' @importFrom gridExtra grid.arrange
#'
#'
#' @examples 
#' # sample 3D p value provided by the package
#' dim(phase2_pval)
#' \donttest{
#' # do the FDR correction
#' pval_fdr = fmri_post_hoc(phase2_pval, 
#'                          fdr_corr = 'fdr',
#'                          spatial_cluster.thr = NULL,
#'                          spatial_cluster.size = NULL, 
#'                          show_comparison = FALSE)
#'
#' # do the spatial clustering
#' pval_posthoc = fmri_post_hoc(pval_fdr,
#'                              fdr_corr = NULL,
#'                              spatial_cluster.thr = 0.05,
#'                              spatial_cluster.size = 5, 
#'                              show_comparison = FALSE)
#' }


fmri_post_hoc = function(p_val_3d,
                         fdr_corr = NULL,
                         spatial_cluster.thr = NULL,
                         spatial_cluster.size = NULL,
                         show_comparison = FALSE, 
                         ...) {
    dim1 = dim(p_val_3d)[1]
    dim2 = dim(p_val_3d)[2]
    dim3 = dim(p_val_3d)[3]
    p_val_3d_original = p_val_3d
    
    gg_list = list()
    
    if (is.null(fdr_corr) == FALSE) {
        dim(p_val_3d) = c(dim1 * dim2 * dim3)
        
        
        p_val_3d = p.adjust(p_val_3d, method = fdr_corr, n = length(p_val_3d))
        dim(p_val_3d) = c(dim1, dim2, dim3)
        
        gg_list[[1]] = fmri_pval_comparison_internal(p_val_3d_original, p_val_3d, names = c("raw p values", "p values after fdr correction"), 
            ...)
    }
    
    
    
    if (is.null(spatial_cluster.size) == FALSE & is.null(spatial_cluster.thr) == FALSE) {
        
        spatial_cluster.filter = cluster.threshold(1 - p_val_3d, level.thr = 1 - spatial_cluster.thr, size.thr = spatial_cluster.size)
        p_val_3d = 1 - ((1 - p_val_3d) * spatial_cluster.filter)
        gg_list[[2]] = fmri_pval_comparison_internal(p_val_3d_original, p_val_3d, names = c("raw p values", "p values after post hoc "), 
            ...)
    }
    
    
    
    if (show_comparison == TRUE) {
        if (length(gg_list) == 2) {
            gridExtra::grid.arrange(gg_list[[1]], gg_list[[2]], nrow = 2)
        } else if (!is.null(gg_list[[1]])) {
            print(gg_list[[1]])
            
        } else if (!is.null(gg_list[[2]])) {
            print(gg_list[[2]])
        }
        
    }
    
    return(p_val_3d)
}


