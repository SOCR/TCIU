## ----warning=FALSE, message=FALSE---------------------------------------------
require(TCIU)

## ----fig.width = 9, fig.align = "center"--------------------------------------
par(mfrow=c(2,2))
image(mask[,30,], asp=1)
image(mask[30,,], asp=1)
image(mask[,,30], asp=1)

## -----------------------------------------------------------------------------
dim(pval1)
summary(pval1)

## -----------------------------------------------------------------------------
fmri_generate = fmri_simulate_func(dim_data = c(64, 64, 40), mask = mask, 
								   ons = c(1, 21, 41, 61, 81, 101, 121, 141), 
								   dur = c(10, 10, 10, 10, 10, 10, 10, 10))

# the output include simulated fMRI data, its mask, 
# the starting time points of the stimulated period and its duration 
# as well as all the stimulated time points
dim(fmri_generate$fmri_data)

## ---- eval=FALSE--------------------------------------------------------------
#  p_simulate_t_test =
#  fmri_stimulus_detect(fmridata= fmri_generate$fmri_data,
#                       mask = fmri_generate$mask,
#                       stimulus_idx = fmri_generate$on_time,
#                       method = "t-test" ,
#                       ons = fmri_generate$ons,
#                       dur = fmri_generate$dur)
#  
#  dim(p_simulate_t_test)
#  summary(p_simulate_t_test)

## -----------------------------------------------------------------------------
# do the FDR correction
pval_fdr = fmri_post_hoc(pval1 , fdr_corr = "fdr",
						 spatial_cluster.thr = NULL,
						 spatial_cluster.size = NULL, 
						 show_comparison = FALSE)

# do the spatial clustering
pval_posthoc = fmri_post_hoc(pval_fdr, fdr_corr = NULL,
							 spatial_cluster.thr = 0.05,
							 spatial_cluster.size = 5, 
							 show_comparison = FALSE)

## ----warning=FALSE------------------------------------------------------------
for(axis in c("x", "y", "z")){
  axis_i = switch(axis, 
                  "x" = {35},
                  "y" = {30},
                  "z" = {22})
  print(fmri_2dvisual(pval2, list(axis, axis_i), 
                      hemody_data=NULL, mask=fmri_generate$mask, 
                      p_threshold = 0.05, legend_show = TRUE, 
                      method = "scale_p",
                      color_pal = "YlOrRd", multi_pranges=TRUE))
}
			

## ----fig.width = 9, fig.align = "center", warning=FALSE-----------------------
fmri_3dvisual(pval2, fmri_generate$mask, 
							p_threshold = 0.05, method="scale_p", multi_pranges=TRUE)$plot

## ----fig.width = 9, fig.align = "center", warning=FALSE-----------------------
# the two p value are the p value generated based on the simulated fMRI
# and the p value saved in the package and finished post hoc test
fmri_pval_comparison_3d(list(pval2, pval_posthoc), mask, 
				                list(0.05, 0.05), list("scale_p", "scale_p"), 
				                multi_pranges=FALSE)


## ----fig.width = 9, warning=FALSE---------------------------------------------
fmri_pval_comparison_2d(list(pval2, pval_posthoc), 
                        list('pval_simulated', 'pval_posthoc'),
                        list(list(35, 33, 22), list(40, 26, 33)), 
                        hemody_data = NULL, 
                        mask = mask, p_threshold = 0.05, 
                        legend_show = FALSE, method = 'scale_p',
                        color_pal = "YlOrRd", multi_pranges=FALSE)

