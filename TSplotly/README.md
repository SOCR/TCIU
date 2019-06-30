# Description of TSplotly package

<a href="http://www.socr.umich.edu/TCIU/"><img align="middle" src="https://github.com/QJoshua/TSplotly/raw/master/TS_temp_page.png"></a>

**Time series data and ARIMA(X) model visualization by plotly [TSplotly]**

Table of contents
=================

<!--ts-->
   * [Table of contents](#table-of-contents)
   * [Overview](#overview)
   * [R package](#r-package)
   * [Installation](#installation)
   * [Vignettes](#vignettes)
   * [Authors](#Authors)
   * [Contact Information](#Contact-information)
<!--te-->


Overview
========

This page is set up to store package, description and related data of an R package **TSplotly**
 TSplotly package provides an effective mechanism to utilize the powerful plotly package for graphing time series data. It contains 4 core functions: 
    TSplot: create plot_ly plot on time series data or fitted ARIMA(X) models.
    ADDline: add lines on existing TSplot objects, as needed.
    GGtoPY: create a convinent way to transform (reformat) ggplot2 datasets into a format that can work on Plot_ly.
    GTSplot: create multiple plot_ly lines (timeseries) based on data frames containing multiple timeseries data

Now this package is mainly used on SOCR's project: ["Time Complexity and Inferential Uncertainty (TCIU)"](http://www.socr.umich.edu/TCIU/). You could see most of the examples in the last part of Chapter 1

R package
=========
The TSplotly protocol has been developed in the [R environment](https://www.r-project.org), see the [TSplotly R package download site on CRAN](https://cran.r-project.org/package=TSplotly) for the latest R version.

To run the TSplotly package. You need to install the following package first: zoo, ggplot2, plotly, forecast.

Installation
============
The version 1.1.0 of the CBDA package can be downloaded and installed with the following command:
```{r Installation of the CBDA package from CRAN, eval = FALSE}
install.packages("TSplotly",repos = 'https://cran.r-project.org/')
```

You could also find the installation package in this repository [here]() or in the Github page of [TCIU project](). 

You could download to your local enviroment and using the following commands to install this package:

```{r Installation of the CBDA package, eval = FALSE}
# Installation from the Windows binary (recommended for Windows systems)
install.packages("~/TSplotly_1.1.0.tar.gz", repos = NULL, type = "win.binary") 
# For windows system, you could also using these following commands:
setwd("~/") ## Set to the working directory containing the "TSplotly_1.1.0.tar.gz" file
system("R CMD INSTALL TSplotly")
# Installation from the source (recommended for Macs and Linux systems)
install.packages("~/TSplotly_1.1.0.tar.gz", repos = NULL, type = "source")
```

Vignettes
=========
The documentation and vignettes, as well as the source and binary files can be found on  [CRAN](https://cran.r-project.org/web/packages/TSplotly/index.html). 

The necessary packages to run the CBDA algortihm are installed automatically at installation.

Authors
===============
This package is created by [Yongkai Qiu](https://socr.umich.edu/people/Yongkai_Qiu.html). Authors of this package includes [Yongkai Qiu](https://socr.umich.edu/people/Yongkai_Qiu.html), [Zhe Yin](http://socr.umich.edu/people/Zhe_Yin.html), and [Ivo D. Dinov](http://www.socr.umich.edu/people/dinov/). This package is set up under the SOCR [TCIU project](http://www.socr.umich.edu/TCIU/). More details of this project would be found on the [SOCR website](http://www.socr.umich.edu/).

Contact Information
===============
[Yongkai Qiu](mailto:yongkai@umich.edu)
