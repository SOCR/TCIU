## TCIU Package

**Spacekime Analytics, Time Complexity and Inferential Uncertainty (TCIU)**

<a href="http://socr.umich.edu/TCIU"><img align="middle" src="https://raw.githubusercontent.com/SOCR/TCIU/master/images/TCUI_P2.png"></a>

### Introduction

The TCIU package provides the core functionality to transform longitudinal data to complex-time (kime) data using analytic and numerical techniques, visualize the original time-series and reconstructed kime-surfaces, perform model based (e.g., tensor-linear regression) and model-free classification and clustering methods.


### Installation

The version 1.0.0 package can be downloaded and installed with the following command once we submitted the package towards CRAN:
```{r Installation of the CBDA package from CRAN, eval = FALSE}
install.packages("TCIU",repos = 'https://cran.r-project.org/')
```

You could also download the [source file](https://github.com/SOCR/TCIU/raw/master/package/TCIU_1.0.0.tar.gz) to your local environment and using the following commands to install this package:

```{r Installation of the CBDA package, eval = FALSE}
# Installation from the binary file
# setwd("...")
# set the working directory to the path where the file is installed
install.packages("./TCIU_1.0.0.tar.gz", repos = NULL, type = "source")
```

### Usage

+ **Validation using 4D Complex (mag,phase) and Real (solely-mag reconstruction) fMRI data**
+ For details, check [TCIU Predictive Analytics](https://www.socr.umich.edu/TCIU/HTMLs/TCIU_Predictive_Analytics.html) for a clear tutorial or refer to the workflow vignette inside the package 
  + **fmri_simulate_func**: an real-valued fMRI data simulation function, used to simply generate a 3 dimension area with activated parts inside
  + **fmri_time_series**: create four interactive time series graphs for the real, imaginary, magnitude and phase parts for the fMRI spacetime data
  + **fmri_kimesurface**: transform the fMRI time-series data at a fixed voxel location into a kimesurface
  + **fmri_image**: create images for front view, side view, and top view of the fMRI image
  + **fmri_ts_forecast**: forecast the fMRI data based on the time series
  + **fmri_stimulus_detect & fmri_post_hoc**: take in real/complex valued fMRI data, detect locations where stimulus occurs and conduct the post-hoc process
  + **fmri_2dvisual & fmri_3dvisual**: a visualization method to display 2d & 3d p-values
  + **fmri_pval_comparison_2d & fmri_pval_comparison_3d**: a visualization method to compare 2d & 3d p-values

+ **Laplace Transform and Kimesurface Transform of TCIU Analytics**
+ For details, check [TCIU Laplace Transform](https://www.socr.umich.edu/TCIU/HTMLs/Laplace_Transform_Timeseries_Kimesurfaces.html) for a clear tutorial or refer to the LT-kimesurface vignettes inside the package 
  + **LT**: numerical method to compute Laplace Transform 
  + **ILT**:  numerical method to compute inverse of Laplace Transform 
  + **kimesurface_transform**: apply the kimesurface transform on a function with a specified set of complex value
  + **inv_kimesurface_transform**:  apply the inverse kimesurface transform to convert to get the original 1D function in [0, 2*pi] or other similar periodic time range


### Contact Information
This package is created by [Yunjie Guo](mailto:jerryguo@umich.edu). Authors of this package includes [Yongkai Qiu](https://socr.umich.edu/people/Yongkai_Qiu.html), [Zhe Yin](http://socr.umich.edu/people/Zhe_Yin.html), Jinwen Cao, Yupeng Zhang, Rongqian Zhang, Yuyao Liu and [Ivo D. Dinov](http://www.socr.umich.edu/people/dinov/). This package is set up under the SOCR [TCIU project](http://www.socr.umich.edu/TCIU/). More details of this project would be found on the [SOCR website](http://www.socr.umich.edu/).

