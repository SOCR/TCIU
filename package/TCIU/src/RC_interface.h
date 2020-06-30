#ifndef RC_INTERFACE_H
#define RC_INTERFACE_H

void read_current_complex(double ***yR_array, double ***yI_array,  
	int i, int j, int n, double *yR, double *yI);

void copy_1d_to_3d(double *vec, double ***array, int d1, int d2, int d3);

void copy_1d_to_2d(double *vec, double **mat, int nrow, int ncol);

void copy_2d_to_1d(double **mat, int nrow, int ncol, double *vec);

void copy_3d_to_1d(double ***array3, int d1, int d2, int d3, double *array1);

void copy_4d_to_1d(double ****array4, int d1, int d2, int d3, int d4, 
	double *array1);

void Rprint_vector(double *a, int length, const char *format);

void Rprint_matrix(double **a, int rows, int cols, const char *format);

void copy_2d_to_1d_int(int **mat, int nrow, int ncol, int *vec);

void copy_1d_to_2d_int(int *vec, int **mat, int nrow, int ncol);

#endif