#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include "array.h"
#include "mat_vec.h"

void read_current_complex(double ***yR_array, double ***yI_array,  
	int i, int j, int n, double *yR, double *yI)
{
	int t;
	for(t=0; t<n; ++t){
		yR[t] = yR_array[t][i][j];
		yI[t] = yI_array[t][i][j];
	}
}

void copy_1d_to_3d(double *vec, double ***array, int d1, int d2, int d3)
{
	int i, j, k, index=0;
	for(i=0; i<d3; ++i){
		for(j=0; j<d2; ++j){
			for(k=0; k<d1; ++k){
				array[k][j][i] = vec[index];
				++index;
			}
		}
	}
}

void copy_1d_to_2d(double *vec, double **mat, int nrow, int ncol)
{
	int i, j, index=0;
	for(j=0; j<ncol;++j){
		for(i=0; i<nrow; ++i){
			mat[i][j] = vec[index];
			++index;
		}
	}
}

void copy_2d_to_1d(double **mat, int nrow, int ncol, double *vec)
{
	int index = 0, i, j;
	for(j=0; j<ncol; ++j){
		for(i=0; i<nrow; ++i){
			vec[index] = mat[i][j];
			++index;
		}
	}
}

void copy_3d_to_1d(double ***array3, int d1, int d2, int d3, double *array1)
{
	int index = 0, i, j, k;
	for(k=0; k<d3; ++k){
		for(j=0; j<d2; ++j){
			for(i=0; i<d1; ++i){
				array1[index] = array3[i][j][k];
				++index;
			}
		}
	}
}

void copy_4d_to_1d(double ****array4, int d1, int d2, int d3, int d4, double *array1)
{
	int index = 0, i, j, k, m;
	for(m=0; m<d4; ++m){
		for(k=0; k<d3; ++k){
			for(j=0; j<d2; ++j){
				for(i=0; i<d1; ++i){
					array1[index] = array4[i][j][k][m];
					++index;
				}
			}
		}
	}
}

void Rprint_vector(double *a, int length, const char *format)
{
	int i;
	for(i=0; i<length; i++)    Rprintf(format, a[i]);
	Rprintf("\n");
}

void Rprint_matrix(double **a, int rows, int cols, const char *format)
{
	int i, j;
	for(i=0; i<rows; i++){
		for(j=0; j<cols; j++) Rprintf(format, a[i][j]);
		Rprintf("\n");
	}
	Rprintf("\n");
}

void copy_2d_to_1d_int(int **mat, int nrow, int ncol, int *vec)
{
	int index = 0, i, j;
	for(j=0; j<ncol; ++j){
		for(i=0; i<nrow; ++i){
			vec[index] = mat[i][j];
			++index;
		}
	}
}

void copy_1d_to_2d_int(int *vec, int **mat, int nrow, int ncol)
{
	int i, j, index=0;
	for(j=0; j<ncol;++j){
		for(i=0; i<nrow; ++i){
			mat[i][j] = vec[index];
			++index;
		}
	}
}