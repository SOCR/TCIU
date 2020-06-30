# @title fmri_on_off_volume
# @description This function calculates volume of On-Off difference period polar.volume of a voxel.
# \strong{Notice:} This function mainly serves as a functional part for other functions. So we
# seldom apply this function by itself.
# 
# @param data a 4D array contains information for the fMRI spacetime image. The data should only contains the magnitude for the fMRI image.
# @param x voxel location for the x axis.
# @param y voxel location for the y axis.
# @param z voxel location for the z axis.
# @param coordinates type of the input position coordinates.
#
# @details The function \code{fmriOnOffVolume.Diff} is used to calculate volume of On-Off difference period polar.volume of a voxel.
# 
# @author SOCR team <\url{http://socr.umich.edu/people/}>
#
# @return the volume of On-Off difference period polar.volume of a voxel.
# @export
# 
# @examples
# fmri_generate = fmri_simulate_func(dim_data = c(64, 64, 40), mask = mask)
# diff = fmriOnOffVolume.Diff(fmri_generate$fmri_data, 40, 41, 33, coordinates = "cartesian")
# 

fmri_on_off_volume = function(data, x, y, z, coordinates = "polar"){
  voxel_data = data[x,y,z,]
  On_data = voxel_data[rep(c(TRUE,FALSE), c(10,10)) ]
  #print(active_data)
  Off_data = voxel_data[rep(c(FALSE,TRUE), c(10,10)) ]
  #print(resting_data)
  
  switch(coordinates, 
         "polar"= {
           volume.diff = 1/80 * ( sum( On_data^2 ) - sum( Off_data^2 ) )
           return(volume.diff)},
         
         "cartesian"= {
           # randomly generate 8 phi kime-phases for each of the 10 time radii
           phi_8_vec <- matrix(NA, ncol=10, nrow = 8)
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
           matrix_ON <- matrix(0, nrow = 21, ncol = 21) 
           matrix_OFF <- matrix(0, nrow = 21, ncol = 21) 
           for (t in 1:10) {
             for (p in 1:8) {
               x = 11+t*cos(phi_8_vec[p,t])
               y = 11+t*sin(phi_8_vec[p,t])
               matrix_ON[x,y]  <- On_data[(p-1)*10 +t]
               matrix_OFF[x,y] <- Off_data[(p-1)*10 +t]
             }
           }
           volume.diff = 1/(21^2) * (sum( matrix_ON^2 ) - sum( matrix_OFF^2 ))
           return( volume.diff )},
         
         {return("Please type a valid coordinate system type!!!")}
  )
  
}