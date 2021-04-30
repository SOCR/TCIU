#' @title interactive graph object of 3D kime-series
#' @description Use \code{plotly} to display in 3D the kime-series as 2D manifolds (kimesurface) over the cartesian domain.
#'
#' @param fmridata a 4d array which contains the spatial and temporal record of fMRI result or a single real valued vector.
#' @param voxel_location a 3d array indicating the spatial location of the brain.
#' @param is.4d The default is true. If change to false, need to input a vector instead of array.
#'
#' @details The function \code{fmri_kimesurface} is display in 3D the kime-series as 2D manifolds (kimesurface) over the Cartesian domain. It helps transform the fMRI time-series data at a fixed voxel location into a kimesurface (kime-series). User can choose to provide the 4D array of the fMRI spacetime image and the voxel_location or a single time-series vector, then a 3D visualization will be shown.
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return an interactive plot in 3D kimesurface
#' @export
#'
#' @import plotly DT scales
#' @importFrom extraDistr rlaplace
#' @importFrom spatstat.core blur
#' @importFrom spatstat.geom as.im
#' 
#' @examples
#' # sample fMRI time-series vector of a single voxel
#' sample_voxel = sample[[5]]
#' fmri_kimesurface(sample_voxel, is.4d = FALSE)[[1]]
#' fmri_kimesurface(sample_voxel, is.4d = FALSE)[[2]]
#' fmri_kimesurface(sample_voxel, is.4d = FALSE)[[3]]
#' fmri_kimesurface(sample_voxel, is.4d = FALSE)[[4]]
#' 

fmri_kimesurface <- function(fmridata, 
                             voxel_location = NULL,
                             is.4d = TRUE) {
  # randomly generate 8 phi kime-phases for each of the 10 time
  phi_8_vec <- matrix(NA, ncol=10, nrow = 8)
  #if(rand_opt=="laplace"){
    for (t in 1:10) { 
      # for a given t, generate 8 new phases
      set.seed(t);
      phi_8_vec[ ,t] <-
        extraDistr::rlaplace(8,mu=0,sigma=0.5)
      # rank-order the phases for consistency
      # within the same foliation leaf
      phi_8_vec[ ,t] <- sort(phi_8_vec[ ,t])
      # force phases in [-pi: pi)
      for (i in 1:8) {
        if (phi_8_vec[i,t] < -pi) 
          phi_8_vec[i,t] <- -pi
        if (phi_8_vec[i,t] >= pi) 
          phi_8_vec[i,t] <- pi
      }
    }
  #}
  if (is.4d == TRUE & is.null(voxel_location) == FALSE ){
    Voxel = fmridata[voxel_location[1],
                     voxel_location[2],
                     voxel_location[3], ]
  }else{
    Voxel = fmridata
  }

  fMRI_ON<-Voxel[c(rep(TRUE,10),rep(FALSE,10))]
  fMRI_OFF<-Voxel[c(rep(FALSE,10),rep(TRUE,10))]
  
  # construct the 160 (time) by 3 (fesatures) DF
  df3D_ON <- data.frame(time=1:10, phi=c(phi_8_vec[1,],phi_8_vec[1,],phi_8_vec[2,],phi_8_vec[2,],phi_8_vec[3,],phi_8_vec[3,],
                                         phi_8_vec[4,],phi_8_vec[4,],phi_8_vec[5,],phi_8_vec[5,],phi_8_vec[6,],phi_8_vec[6,],
                                         phi_8_vec[7,],phi_8_vec[7,],phi_8_vec[8,],phi_8_vec[8,]), switch=c(rep(TRUE,10),rep(FALSE,10)), fMRI=Voxel)
  # dim(df3D_ON); head(df3D_ON, 15)
  
  # Convert the long DF representing fMRI_ON and fMRI_OFF from polar coordinates to Cartesian coordinates
  matrix_ON <- matrix(0, nrow = 21, ncol = 21) 
  matrix_OFF <- matrix(0, nrow = 21, ncol = 21) 
  for (t in 1:10) {
    for (p in 1:8) {
      x = 11+t*cos(phi_8_vec[p,t])
      y = 11+t*sin(phi_8_vec[p,t])
      matrix_ON[x,y]  <- fMRI_ON[(p-1)*10 +t]
      matrix_OFF[x,y] <- fMRI_OFF[(p-1)*10 +t]
    }
  }
  
  # fix the plot_ly Text Lables
  x <- vector()
  y <- vector()
  i <- 1
  for (t in 1:10) {
    for (p in 1:8) {
      x[i] = 11+t*cos(phi_8_vec[p,t])
      y[i] = 11+t*sin(phi_8_vec[p,t])
      i <- i+1
    }
  }
  
  matrix_ON_smooth <- (1/10000)*as.matrix(spatstat.core::blur(spatstat.geom::as.im(matrix_ON), sigma=0.5))
  matrix_OFF_smooth <- (1/10000)*as.matrix(spatstat.core::blur(spatstat.geom::as.im(matrix_OFF), sigma=0.5))
  
  hoverText <- cbind(x=1:21, y=1:21, height=as.vector(t(matrix_ON_smooth))) # tail(mytext)
  custom_txt <- matrix(NA, nrow=21, ncol=21)
  hoverTextOFF <- cbind(x=1:21, y=1:21, height=as.vector(t(matrix_OFF_smooth))) # tail(mytext)
  custom_txtOFF <- matrix(NA, nrow=21, ncol=21)
  
  for (x in 1:21) {
    for (y in 1:21) {
      t = sqrt((x-11)^2 + (y-11)^2)
      p = atan2(y-11, x-11)
      custom_txt[x,y] <- paste(' fMRI: ', round(hoverText[(x-1)*21+y, 3], 3),
                               '\n time: ', round(t, 0),
                               '\n phi: ', round(p, 2)) 
      custom_txtOFF[x,y] <- paste(' fMRI: ', round(hoverTextOFF[(x-1)*21+y, 3], 3),
                                  '\n time: ', round(t, 0),
                                  '\n phi: ', round(p, 2)) 
    }
  }
  
  xx2 <- 11 + c(-10:10) %o% cos(seq(-pi, pi, 2*pi/20))
  yy2 <- 11 + c(-10:10) %o% sin(seq(-pi, pi, 2*pi/20))
  #zz2 <- as.vector(matrix_ON_smooth) %o% rep(1, 21*21)
  zz2 <- matrix_ON_smooth
  ww2 <- matrix_OFF_smooth
  dd2 <- matrix_ON_smooth-matrix_OFF_smooth
  
  dd2scale<-fmri_split_ab_bl(dd2)
  
  #plot 2D into 3D and make the text of the diameter (time), height (r), and phase (phi)
  f <- list(family = "Courier New, monospace", size = 18, color = "black")
  x <- list(title = "k1", titlefont = f)
  y <- list(title = "k2", titlefont = f)
  z <- list(title = "fMRI Kime-series", titlefont = f)
  zd <- list(title = "fMRI Kime-ON/OFF difference", titlefont = f)
  
  plot1<-plot_ly(x = ~xx2, y = ~yy2, z = ~zz2, type = "surface", colors=c("#FFFFFF","#0000FF"), # scatterpolar
          text = custom_txt, hoverinfo = "text", showlegend = FALSE) %>% 
    #add_trace(x=~xx2, y=~yy2, z=~ww2, colors = c("blue", "yellow"),
    #          type="surface", text = custom_txtOFF, hoverinfo = "text",
    #          opacity=0.3, showscale = FALSE, showlegend = FALSE) %>%
    # trace the main Z-axis
    add_trace(x=11, y=11, z=0:0.15, type="scatter3d", mode="lines", 
              line = list(width = 10, color="red"), name="Space(x)", 
              hoverinfo="none", showlegend = FALSE) %>%
    layout(dragmode = "turntable", title = "ON Kime-Surface/Kime-Series at a fixed voxel location",
           scene = list(xaxis = x, yaxis = y, zaxis = z), showlegend = FALSE)
  
  # Plot OFF kime-surface
  plot2<-plot_ly(x = ~xx2, y = ~yy2, z = ~ww2, type = "surface",colors=c("#FFFFFF","#0000FF"),   # scatterpolar
          text = custom_txt, hoverinfo = "text", showlegend = FALSE) %>% 
    #add_trace(x=~xx2, y=~yy2, z=~ww2, colors = c("blue", "yellow"),
    #          type="surface", text = custom_txtOFF, hoverinfo = "text",
    #          opacity=0.3, showscale = FALSE, showlegend = FALSE) %>%
    # trace the main Z-axis
    add_trace(x=11, y=11, z=0:0.15, type="scatter3d", mode="lines", 
              line = list(width = 10, color="red"), name="Space(x)", 
              hoverinfo="none", showlegend = FALSE) %>%
    layout(dragmode = "turntable", title = "OFF Kime-Surface/Kime-Series at a fixed voxel location",
           scene = list(xaxis = x, yaxis = y, zaxis = z), showlegend = FALSE)
  
  plot3<-plot_ly(x = ~xx2, y = ~yy2, z = ~dd2, type = "surface",colors = dd2scale, #colors=c("#FFFF00","#FFFFFF","#0000FF"),  # scatterpolar
          text = custom_txt, hoverinfo = "text", showlegend = FALSE) %>% 
    #add_trace(x=~xx2, y=~yy2, z=~ww2, colors = c("blue", "yellow"),
    #          type="surface", text = custom_txtOFF, hoverinfo = "text",
    #          opacity=0.3, showscale = FALSE, showlegend = FALSE) %>%
    # trace the main Z-axis
    add_trace(x=11, y=11, z=-0.15:0.15, type="scatter3d", mode="lines", 
              line = list(width = 10, color="red"), name="Space(x)", 
              hoverinfo="none", showlegend = FALSE) %>%
    layout(dragmode = "turntable", title = "Difference for ON & OFF Kime-Surface/Kime-Series at a fixed voxel location",
           scene = list(xaxis = x, yaxis = y, zaxis = zd), showlegend = FALSE)
  
  return(list(df3D_ON,plot1,plot2,plot3))
}