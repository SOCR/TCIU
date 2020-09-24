# TCIU

**Data Science, Time Complexity and Inferential Uncertainty (TCIU)**

<a href="http://socr.umich.edu/TCIU"><img align="middle" src="https://raw.githubusercontent.com/SOCR/TCIU/master/images/TCUI_P1.png"></a>

Table of contents
=================

<!--ts-->
   * [Table of contents](#table-of-contents)
   * [Overview](#overview)
   * [R Code](#r-code)
   * [Team](#team)
   * [Acknowledgments](#acknowledgments)
   * [References](#references)
<!--te-->


Overview
========

The [SOCR Data Science Fundamentals project explores new theoretical representation and analytical strategies to understand large and complex data, time complexity and inferential uncertainty](http://www.socr.umich.edu/TCIU/). It utilizes information measures, entropy KL divergence, PDEs, Dirac’s bra-ket operators (&#12296; , &#12297;). This fundamentals of data science research project will explore time-complexity and inferential uncertainty in modeling, analysis and interpretation of large, heterogeneous, multi-source, multi-scale, incomplete, incongruent, and longitudinal data. 

R Code
======

The examples, demonstrations and simulations are designed, built, implemented and validated in the [R environment](https://www.r-project.org). 

The source R code for the package is in the [package folder](https://github.com/SOCR/TCIU/tree/master/package).

The source RMarkdown code for the website is in the [website folder](https://github.com/SOCR/TCIU/tree/master/website).

To interactively run of the demo code, following a [TCIU package installation](https://github.com/SOCR/TCIU/tree/master/package/TCIU), use the following command in the R console/shell:
```{r}
library(TCIU)
demo(fmri_demo_func, package="TCIU")
```

Team
====

[SOCR Team](http://www.socr.umich.edu/people/) including [Ivo D. Dinov](http://umich.edu/~dinov), [Milen V. Velev]( http://www.socr.umich.edu/people/Milen_Velev.html), [Yongkai Qiu]( http://www.socr.umich.edu/people/Yongkai_Qiu.html), [Zhe Yin]( http://www.socr.umich.edu/people/Zhe_Yin.html), Yufei Yang, Yunjie Guo, Yupeng Zhang, Rongqian Zhang, Yuyao Liu, Jinwen Cao, Zijing Li, Daxuan Deng, Yueyang Shen, and others.

Acknowledgments
===============

This work is supported in part by NIH grants [P20 NR015331](www.socr.umich.edu/CSCD), [UL1TR002240](https://projectreporter.nih.gov/project_info_description.cfm?aid=9491961&icde=39078316), [P30 DK089503](http://mmoc.med.umich.edu/), [UL1TR002240](https://www.michr.umich.edu), and NSF grants [1916425](http://midas.umich.edu/), [1734853](http://brain-life.org/), [1636840](http://neurosciencenetwork.org/), [1416953](http://distributome.org), [0716055](http://socr.umich.edu) and [1023115](http://distributome.org). Students, trainees, scholars, and researchers from SOCR, BDDS, MNORC, MIDAS, MADC, MICHR, and the broad R-statistical computing community have contributed ideas, code, and support.

References
==========

* Dinov, ID and Velev, MV (2021) [Data Science: Time Complexity, Inferential Uncertainty, and Spacekime Analytics](https://www.degruyter.com/view/title/576646), De Gruyter STEM Series, Berlin/Boston, ISBN 9783110697803 / 3110697807. 
* [TCIU Website](https://tciu.predictive.space/).
* [Spacekime](https://Spacekime.org).
* [TSPlotly package on CRAN](https://cran.r-project.org/web/packages/TSplotly/index.html).
