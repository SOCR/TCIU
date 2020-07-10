## ----warning=FALSE, message=FALSE---------------------------------------------
require(TCIU)
require(doParallel)
require(cubature)
require(oro.nifti)
require(magrittr)
require(plotly)
require(ggplot2)

## ----eval=FALSE, warning=FALSE, message=FALSE, fig.align = "center"-----------
#  # For this part of code, we comment it out and import the output plot already generated before to reduce time
#  # But this part of code can be run successfully. If you are interested, you can try it on your computer!
#  
#  # discrete Laplace Transform of sine
#  range_limit = 2
#  x2 = seq(from = 0, to = range_limit, length.out = 50)[2:50]
#    # drop the first row to avoid real part value of 0
#  y2 = seq(from = 0, to = range_limit, length.out = 50)[2:50]
#    # drop the first column to avoid imaginary part value of 0
#  
#  # Recompute the LT(sin) discretized on lower-res grid
#  z2_grid = array(dim=c(length(x2), length(y2)))# x2 %o% y2
#  
#  f_sin <- function(t) { sin(t) }
#  
#  # kime surface transform
#  # use parallel computing to speed up code
#  ncors = 2 # please choose the ncors according to the number of cores your PC has
#  # it is better that you increase the number of cores used for parallel computing if your computer allows
#  cl <- makeCluster(ncors)
#  registerDoParallel(cl)
#  F = list()
#  for (i in 1:length(x2) ){
#    F[[i]] =
#      foreach(j = 1:length(y2),
#              .export='cubintegrate',
#              .packages='cubature') %dopar% {
#        TCIU::LT(FUNCT=f_sin, complex(real=x2[i], imaginary = y2[j]))
#      }
#  }
#  
#  stopCluster(cl)
#  F_vec = lapply(F, unlist)
#  z2_grid = unlist(do.call(rbind, F_vec))
#  
#  # explicit form of Laplace Transform of sine
#  laplace_sine = function(p) { 1/(p^2 + 1) } # Exact laplace transform of sin(x), continuous function
#  
#  XY = expand.grid(X=x2, Y=y2)
#  complex_xy = mapply(complex, real=XY$X,imaginary=XY$Y)
#  sine_z =laplace_sine(complex_xy)
#  dim(sine_z) = c(length(x2), length(y2))# dim(sine_z)  # [1] 49 49
#  
#  # make the two plots in the same plot to compare
#  lt_ilt_plotly = plot_ly(hoverinfo="none", showscale = FALSE)%>%
#    add_trace(z=Re(sine_z)-1, type="surface", surfacecolor=Im(sine_z))  %>%
#    add_trace(z = Re(z2_grid), type="surface", opacity=0.7, surfacecolor=Im(z2_grid) )%>%
#    layout(title =
#    "Laplace Transform, LT(sin()), Height=Re(LT(sin())), Color=Re(LT(sin())) \n Contrast Exact (Continuous) vs.
#           Approximate (Discrete) Laplace Transform", showlegend = FALSE)
#  
#  lt_ilt_plotly
#  

## ----fig.align = "center", warning=FALSE, message=FALSE-----------------------
lt_ilt_plotly = system.file("extdata", "lt_ilt_plotly.rda", package = "TCIU", mustWork=TRUE)
get(load(lt_ilt_plotly))

## ----warning=FALSE, message=FALSE, fig.width = 9, fig.align = "center", eval=FALSE----
#  # For this part of code, we comment it out and import the output plot already generated before to reduce time
#  # But this part of code can be run successfully. If you are interested, you can try it on your computer!
#  
#  # discrete Laplace Transform of sine
#  f_sin = function(t) { sin(t) }
#  lt_sine = function(z) TCIU::LT(f_sin, z)
#  
#  # inverse Laplace Transform on the lt_sine
#  # using parallel computing to speed up code
#  tvalsn <- seq(0, pi*2, length.out = 20)
#  cl <- makeCluster(ncors)
#  registerDoParallel(cl)
#  sinvalsn <- foreach(t=1:length(tvalsn),
#                      .export='cubintegrate',
#                      .packages='cubature')  %dopar% {
#    TCIU::ILT(FUNCT=lt_sine, t=tvalsn[t])
#    }
#  stopCluster(cl)
#  sinvalsn = unlist(sinvalsn)
#  
#  # make the plot of the result from ILT
#  # to see whether it still looks like sine
#  sinvalsn_df2 <- as.data.frame(cbind(Re=Re(sinvalsn),Im=Im(sinvalsn),
#                                      Sin=sin(tvalsn), time_points=tvalsn))
#  lt_ilt_sine = ggplot(sinvalsn_df2, aes(x=time_points))+
#    geom_line(aes(y=Re, color="Real"), linetype=1, lwd=2) +
#    geom_line(aes(y = Sin, color="Sin"), linetype=2, lwd=1) +
#    scale_color_manual(name="Index",
#                       values = c("Real"="steelblue", "Sin"="darkred"))+
#    labs(title = "Original fMRI Time-series f(t)=sin(t) and \n Reconstructed f'(t)=ILT(F)=ILT(discrete LT(f))",
#         subtitle = bquote("F" ~ "=" ~ "discrete LT(sine)")) +
#    xlab("Time") + ylab("fMRI Image Intensities (f and f')") +
#    theme_grey(base_size = 16) +
#    theme(legend.title = element_text(size=14, color = "black", face="bold"),
#          panel.grid.minor.y = element_blank(),
#          panel.grid.major.y = element_blank(),
#          plot.title = element_text(hjust = 0.5),
#          plot.subtitle = element_text(hjust = 0.5))
#  
#  lt_ilt_sine

## ---- fig.align = "center", warning=FALSE, message=FALSE, fig.width = 9-------
lt_ilt_sine = system.file("extdata", "lt_ilt_sine.rda", package = "TCIU", mustWork=TRUE)
get(load(lt_ilt_sine))

## ----warning=FALSE, message=FALSE---------------------------------------------
x = seq(0, 2, length.out=50)[2:50]; y = seq(0, 2, length.out=50)[2:50];
# do Kimesurface transform on sine function
z2_grid = kimesurface_transform(FUNCT = function(t) { sin(t) },
                      real_x = x, img_y = y)

# make the plot after Kimesurface transformation
surf_color <- atan2(Im(z2_grid), Re(z2_grid))
colorscale = cbind(seq(0, 1, by=1/(length(x) - 1)), rainbow(length(x)))
magnitude <- (sqrt( Re(z2_grid)^2+ Im(z2_grid)^2))


p <- plot_ly(hoverinfo="none", showscale = FALSE) %>%
    add_trace(z = magnitude, 
              surfacecolor=surf_color, colorscale=colorscale,   #Phase-based color
              type = 'surface', opacity=1, visible=T) %>%
    layout(title = "fMRI Kime-Surface, F=LT(fMRI) \n Height=Mag(F), Color=Phase(F)", showlegend = FALSE,
           scene = list(aspectmode = "manual", aspectratio = list(x=1, y=1, z=1.0) ) ) # 1:1:1 aspect ratio
p

## ----warning=FALSE, message=FALSE, eval = FALSE-------------------------------
#  # load the fMRI data
#  
#  fMRIURL <- "http://socr.umich.edu/HTML5/BrainViewer/data/fMRI_FilteredData_4D.nii.gz"
#  fMRIFile <- file.path(tempdir(), "fMRI_FilteredData_4D.nii.gz")
#  download.file(fMRIURL, dest=fMRIFile, quiet=TRUE)
#  fMRIVolume <- readNIfTI(fMRIFile, reorient=FALSE)
#  # dimensions: 64 x 64 x 21 x 180 ; 4mm x 4mm x 6mm x 3 sec
#  
#  
#  # extract a set of time series data from a voxel of fMRI data
#  xA_fMRI_1D_x20_y20_z11 <- fMRIVolume[20, 20, 11, ]
#  
#  # get the smooth function of this time series data
#  # Instead of using the extremely noisy fMRI data and avoiding integration problems,
#  # smooth "f" and use the **smooth version, f**
#  time_points <- seq(0+0.001, 2*pi, length.out = 180)
#  f <- smooth.spline(ceiling((180*time_points)/(2*pi)),
#                     xA_fMRI_1D_x20_y20_z11, df = 10)$y
#  
#  # Define the f(t)=smooth(fMRI)(t) signal as a function of real time 0<t<=2*pi
#  fmri_funct <- function(t) {
#    if (t < 0+0.001 || t > 2*pi) {  return ( 0 )  } else {
#      return ( f[ceiling((180*t)/(2*pi))] )
#    }
#  }
#  
#  
#  # do the Kimesurface transform
#  ncors = 8
#  # please choose the ncors according to the number of cores your PC has
#  x = seq(0, 2, length.out=50)[2:50]; y = seq(0, 2, length.out=50)[2:50];
#  # do Kimesurface transform on sine function
#  z2_grid_fmri = kimesurface_transform(FUNCT = fmri_funct,
#                                       glb_para="f",
#                                       real_x = x, img_y = y,
#                                       parallel_computing = TRUE,
#                                       ncor=ncors)
#  
#  # make the plot of function after Kimesurface transformation
#  surf_color <- atan2(Im(z2_grid_fmri), Re(z2_grid_fmri))
#  colorscale = cbind(seq(0, 1, by=1/(length(x) - 1)), rainbow(length(x)))
#  magnitude <- (sqrt( Re(z2_grid_fmri)^2+ Im(z2_grid_fmri)^2))
#  
#  
#  p_fmri <- plot_ly(hoverinfo="none", showscale = FALSE) %>%
#      add_trace(z = magnitude,
#                surfacecolor=surf_color, colorscale=colorscale,   #Phase-based color
#                type = 'surface', opacity=1, visible=T) %>%
#      layout(title = "fMRI Kime-Surface, F=LT(fMRI) \n Height=Mag(F), Color=Phase(F)", showlegend = FALSE,
#             scene = list(aspectmode = "manual", aspectratio = list(x=1, y=1, z=1.0) ) ) # 1:1:1 aspect ratio
#  p_fmri

## ---- fig.align = "center", warning=FALSE, message=FALSE, fig.width = 9-------
p_fmri = system.file("extdata", "p_fmri.rda", package = "TCIU", mustWork=TRUE)
get(load(p_fmri))

## ----warning=FALSE, message=FALSE, fig.width = 9, fig.align = "center"--------
time_points <- seq(0+0.001, 2*pi, length.out = 180)
inv_data = inv_kimesurface_transform(time_points, z2_grid)
inv_data = inv_kimesurface_transform(time_points, z2_grid,num_length = 23,
                                     m=1, msg=TRUE)

time_Intensities_ILT_df2 <- as.data.frame(cbind(Re=scale(Re(inv_data$Smooth_Reconstruction)),
                                                Im=scale(Re(inv_data$Raw_Reconstruction)),
                                                fMRI=scale(Re(sin(time_points))),
                                                time_points=time_points))
colnames(time_Intensities_ILT_df2) = c("Smooth Reconstruction",
                                       "Raw Reconstruction",
                                       "Original Sin", "time_points")
df = reshape2::melt(time_Intensities_ILT_df2, id.var = "time_points")
ggplot(df, aes(x = time_points, y = value, colour = variable)) +
  geom_line(linetype=1, lwd=3) +
  ylab("Function Intensity") + xlab("Time") +
  theme(legend.position="top")+
  labs(title= bquote("Comparison between" ~ "f(t)=Smooth(Sin)(t)" ~ "and Smooth(ILT(LT(Sin)))(t); Range [" ~ 0 ~":"~ 2*pi~"]"))

