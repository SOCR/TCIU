#ifndef MAT_VEC_H
#define MAT_VEC_H

void copy_vec(int len, double *a, double *b);

void outer_prod(int lena, int lenb, double *a, double *b, 
	double **ab);

void make_identity_mat(int n, double **I_n);

double dot_prod(int len, double *a, double *b);
int max_vec(double *vec, int n);
int min(int a, int b);
int max(int a, int b);
void zero_mat(int nrow, int ncol, double **mat);
void OLS(double **X, double *y, int n, int q, double *beta);
double samp_var(double **X, double *r, int n, int q, double *beta);
double dist_max(double *a, double *b, int n);
void par_update(int p, int q, double *alpha, double *beta, double sig2, double *par);
void X_transpose_X(double **X, int nrow, int ncol, double **XpX);
void X_transpose_y(double **X, double *y, int nrow, int ncol, double *Xpy);
int multiply(double **a, int arows, int acols,
	double **b, int brows, int bcols, double **c);
int matxvec(double **a, int arows, int acols,
	double *x, int xrows, double *y);
void transpose(int nrow, int ncol, double **M, double **tM);
void print_dmatrix(double **a, int rows, int cols, const char *format);
void print_dvector(double *a, int rows, const char *format);
double quadratic(double **A, double *x, int p); 
     /*calculates x'Ax where A is a p x p-dimensional matrix */
int cpy(double **a,int nrows,int ncols,double **b);
double bilinear(double *a, double **M, double *b, int n);
	/* Computes computes the bilinear form t(a)%*%M%*%b */
void my_inv(int size, double **mat);
double my_det(int size, double **mat);
int matinv(int sizeA,double **A,double (*determinant));



#endif /* MATVEC_H */
