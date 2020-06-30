#ifndef SPATIAL_SMOOTH_H
#define SPATIAL_SMOOTH_H

void Rwrapper_spatial_smooth3d(double *slice_in_vec, 
	int *n1, int *n2, int *n3, double *kern_vec, int *kerndim);

void spatial_smooth3d(double ***slice_in, int n1, int n2, 
	int n3, double ***kern, int kerndim);

void spatial_smooth1pt_3d(int i, int j, int k, int n1, int n2, 
	int n3, double ***slice_in, double ***kern, int kerndim, 
	int kw, double *smoothed_pt);

void Rwrapper_spatial_smooth2d(double *slice_in_vec, 
	int *n1, int *n2, double *kern_vec, int *kerndim);

void spatial_smooth2d(double **slice_in, int n1, int n2, 
	double **kern, int kerndim);

void spatial_smooth1pt(int i, int j, int n1, int n2, double **slice_in, 
	double **kern, int kerndim, int kw, double *smoothed_pt);

#endif