#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "array.h"
#include "mat_vec.h"
#include "RC_interface.h"
#include "spatial_smooth.h"

void Rwrapper_isnan(int *n, double *x, int *well){
	int i;
	for(i=0; i<*n; ++i){
		if(!isnan(x[i])) well[i] = 1;
		else well[i] = 0;
	}
}

void Rwrapper_spatial_smooth3d(double *slice_in_vec, 
	int *n1, int *n2, int *n3, double *kern_vec, int *kerndim) 
	//NAOK = T
{
	double ***slice_in;	MAKE_3ARRAY(slice_in, *n1, *n2, *n3);
	double ***kern;		MAKE_3ARRAY(kern, *kerndim, *kerndim, *kerndim);
	copy_1d_to_3d(slice_in_vec, slice_in, *n1, *n2, *n3);
	copy_1d_to_3d(kern_vec, kern, *kerndim, *kerndim, *kerndim);
	spatial_smooth3d(slice_in, *n1, *n2, *n3, kern, *kerndim);
	copy_3d_to_1d(slice_in, *n1, *n2, *n3, slice_in_vec);
	FREE_3ARRAY(slice_in);
	FREE_3ARRAY(kern);
}

void spatial_smooth3d(double ***slice_in, int n1, int n2, 
	int n3, double ***kern, int kerndim)
{
	int i, j, k, kw = (kerndim - 1) / 2;
	double ***slice_out;	MAKE_3ARRAY(slice_out, n1, n2, n3);
	for(i=0; i<n1; ++i){
		for(j=0; j<n2; ++j){
			for(k=0; k<n3; ++k){
				if(!isnan(slice_in[i][j][k]))
					spatial_smooth1pt_3d(i, j, k, n1, n2, n3, 
						slice_in, kern, kerndim, kw, 
						&slice_out[i][j][k]);
			}
		}
	}
	for(i=0; i<n1; ++i)
		for(j=0; j<n2; ++j)
			for(k=0; k<n3; ++k)
				if(!isnan(slice_in[i][j][k]))
					slice_in[i][j][k] = slice_out[i][j][k];
	FREE_3ARRAY(slice_out);
}

void spatial_smooth1pt_3d(int i, int j, int k, int n1, int n2, 
	int n3, double ***slice_in, double ***kern, int kerndim, 
	int kw, double *smoothed_pt)
{
	int l, m, p; //indices of kernel
	int inrange, inbrain;
	double sum_kern = 0.0, sum_image = 0.0;
	for(l=0; l<kerndim; ++l){
		for(m=0; m<kerndim; ++m){
			for(p=0; p<kerndim; ++p){
				inrange = ((i+l-kw)>=0) && ((i+l-kw)<n1) && 
					((j+m-kw)>=0) && ((j+m-kw)<n2) &&
					((k+p-kw)>=0) && ((k+p-kw)<n3);
				if(inrange){
					inbrain = !isnan(slice_in[i+l-kw][j+m-kw][k+p-kw]);
					if(inbrain){
						sum_image += slice_in[i+l-kw][j+m-kw][k+p-kw] * 
										kern[l][m][p];
						sum_kern += kern[l][m][p];
					}
				}
			}
		}
	}
	*smoothed_pt = sum_image / sum_kern;
}

void Rwrapper_spatial_smooth2d(double *slice_in_vec, 
	int *n1, int *n2, double *kern_vec, int *kerndim) 
	//NAOK = T
{
	double **slice_in;	MAKE_MATRIX(slice_in, *n1, *n2);
	double **kern;		MAKE_MATRIX(kern, *kerndim, *kerndim);
	copy_1d_to_2d(slice_in_vec, slice_in, *n1, *n2);
	copy_1d_to_2d(kern_vec, kern, *kerndim, *kerndim);
	spatial_smooth2d(slice_in, *n1, *n2, kern, *kerndim);
	copy_2d_to_1d(slice_in, *n1, *n2, slice_in_vec);
	FREE_MATRIX(slice_in);
	FREE_MATRIX(kern);
}

void spatial_smooth2d(double **slice_in, int n1, int n2, 
	double **kern, int kerndim)
{
	int i, j, kw = (kerndim - 1) / 2;
	double **slice_out;		MAKE_MATRIX(slice_out, n1, n2);
	for(i=0; i<n1; ++i){
		for(j=0; j<n2; ++j){
			if(!isnan(slice_in[i][j])) 
				spatial_smooth1pt(i, j, n1, n2, slice_in, kern, 
					kerndim, kw, &slice_out[i][j]);
		}
	}
	for(i=0; i<n1; ++i)
		for(j=0; j<n2; ++j)
			if(!isnan(slice_in[i][j]))
				slice_in[i][j] = slice_out[i][j];
	FREE_MATRIX(slice_out);
}

void spatial_smooth1pt(int i, int j, int n1, int n2, double **slice_in, 
	double **kern, int kerndim, int kw, double *smoothed_pt)
{
	int k, m;
	double sum_kern = 0.0, sum_slice = 0.0;
	for(k=0; k<kerndim; ++k){
		for(m=0; m<kerndim; ++m){
			if(((i+k-kw)>=0) & ((i+k-kw)<n1)
			& ((j+m-kw)>=0) & ((j+m-kw)<n2)){
				if(!isnan(slice_in[i+k-kw][j+m-kw])){
					sum_slice += slice_in[i+k-kw][j+m-kw] * kern[k][m];
					sum_kern += kern[k][m];
				}
			}
		}
	}
	*smoothed_pt = sum_slice / sum_kern;
}
