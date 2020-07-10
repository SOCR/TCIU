## ----warning=FALSE, message=FALSE---------------------------------------------
require(TCIU) 
require(DT)
require(AnalyzeFMRI)

## ----warning=FALSE, message=FALSE---------------------------------------------
# simulate the fMRI data
fmri_generate = fmri_simulate_func(dim_data = c(64, 64, 40), mask = mask, 
                                   ons = c(1, 21, 41, 61, 81, 101, 121, 141), 
                                   dur = c(10, 10, 10, 10, 10, 10, 10, 10))
dim(fmri_generate$fmri_data)

## ----warning=FALSE, message=FALSE---------------------------------------------
sample_voxel = system.file("extdata", "sample_voxel.rda", package = "TCIU", mustWork=TRUE)
reference_plot = system.file("extdata", "reference_plot.rda", package = "TCIU", mustWork=TRUE)
load(sample_voxel)
load(reference_plot)
fmri_time_series(sample_voxel, voxel_location = NULL, is.4d = FALSE, ref = reference_plot)

## ----warning=FALSE, message=FALSE---------------------------------------------
 # a data-frame with 160 rows and 4 columns: time (1:10), phases (8), states (2), and fMRI data (Complex or Real intensity)
 datatable(fmri_kimesurface(fmri_generate$fmri_data, c(44,30,33))[[1]])
 # ON Kime-Surface
 fmri_kimesurface(fmri_generate$fmri_data, c(44,30,33))[[2]]
 # OFF Kime-Surface 
 fmri_kimesurface(fmri_generate$fmri_data, c(44,30,33))[[3]]
 # ON&OFF Kime-Surface 
 fmri_kimesurface(fmri_generate$fmri_data, c(44,30,33))[[4]]

## ----warning=FALSE, message=FALSE---------------------------------------------
fmri_image(fmri_generate$fmri_data, option="manually", voxel_location = c(40,22,33), time=4)

## ----warning=FALSE, message=FALSE---------------------------------------------
smoothmod<-GaussSmoothArray(fmri_generate$fmri_data, sigma = diag(3,3))
fmri_ts_forecast(smoothmod, voxel_location=c(41,44,33))

## ---- eval=FALSE--------------------------------------------------------------
#  # If the dimension of fmridata is lower than 4, the corresponding dimension of output p-values will be lower.
#  # Users can also choose other test method or apply FDR correction and spatial clustering to get the result.
#  # Here, we provide two examples of t-test and wilcoxon-test for the 4D real-valued fMRI data.
#  p_simulate_t_test = fmri_stimulus_detect(fmridata= fmri_generate$fmri_data,
#                                           mask = fmri_generate$mask,
#                                           stimulus_idx = fmri_generate$on_time,
#                                           method = "t-test")
#  p_simulate_wilcoxon = fmri_stimulus_detect(fmridata= fmri_generate$fmri_data,
#                                           mask = fmri_generate$mask,
#                                           stimulus_idx = fmri_generate$on_time,
#                                           method = "wilcoxon-test")

## -----------------------------------------------------------------------------
# load 3D p value stored in the package
dim(pval1)
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

## -----------------------------------------------------------------------------
fmri_2dvisual(pval1, list('x',40), hemody_data=NULL, mask=mask, p_threshold=0.05)
fmri_2dvisual(pval1, list('y',23), hemody_data=NULL, mask=mask, p_threshold=0.05)
fmri_2dvisual(pval1, list('z',33), hemody_data=NULL, mask=mask, p_threshold=0.05)

## -----------------------------------------------------------------------------
fmri_pval_comparison_2d(list(pval1, pval2), 
                        list('pval', 'pval2'),
                        list(list(40, 26, 33), list(40, 26, 33)), 
                        hemody_data = NULL, 
                        mask = mask, p_threshold = 0.05, 
                        legend_show = FALSE, method = 'scale_p',
                        color_pal = "YlOrRd", multi_pranges=TRUE)    

## -----------------------------------------------------------------------------
fmri_3dvisual(pval1, mask, p_threshold = 0.05, method="scale_p")$plot

## -----------------------------------------------------------------------------
fmri_pval_comparison_3d(list(pval1, pval2), mask, 
                        list(0.05, 0.05), list('scale_p', 'scale_p'), multi_pranges=FALSE)

