/* mat_vec.c
 * 
 * Basic vector and matrix functions
 * 
 * Authors: Rouben Rostamian<rostamian@umbc.edu> and Ranjan Maitra<maitra@iastate.edu>
 * Fall 1996
 * Revised January 1999
 * Revised October 2000
 * Thoroughly revised March 2005
*/

#include <stdio.h>		/* for fprintf() */
#include <stdlib.h>		/* for malloc() */
#include <math.h>		/* for sqrt() */
#include <assert.h>
#include "array.h"
#include "mat_vec.h"

void copy_vec(int len, double *a, double *b)
{
	int i;
	for(i=0; i<len; ++i) b[i] = a[i];
}

void outer_prod(int lena, int lenb, double *a, double *b, 
	double **ab)
// Computes a %*% t(b)
{
	int i, j;
	for(i=0; i<lena; ++i)
		for(j=0; j<lenb; ++j)
			ab[i][j] = a[i] * b[j];
}

void make_identity_mat(int n, double **I_n)
{
	int i, j;
	for(i=0; i<n; ++i){
		for(j=0; j<n; ++j){
			if(i==j)	I_n[i][j] = 1.0;
			else	I_n[i][j] = 0.0;
		}
	}
}

double dot_prod(int len, double *a, double *b)
{
	int i;
	double temp = 0.0;
	for(i=0; i<len; ++i)
		temp += a[i] * b[i];
	return temp;
}
	

void Rwrapper_generate_normal(int *n, double *mu, double *sigma, 
	double *x)
{
	int i;
	for(i=0; i<*n; ++i){
		x[i] = generate_normal(*mu, *sigma);
	}
}

void arp_sim(int n, double sig2, double *alpha, double *y, int p)
{
	int init = 100;
	int i, j;
	double temp = 0.0;
	double sigma = sqrt(sig2);
	double *z; 		MAKE_VECTOR(z, n+init);
	double *x;		MAKE_VECTOR(x, n+init);
	for(i=0; i< n+init; i++)		z[i] = generate_normal(0., sigma);			
	for(i=0; i<p; i++)			x[i] = z[i];
	for(i=p; i<n+init; i++){
		for(j=1; j<=p; j++)		temp += alpha[j-1] * x[i-j];
		x[i] = temp + z[i];
		temp= 0.0;
	}
	for(i=0; i<n; i++)		y[i] = x[init+i];
	FREE_VECTOR(z);
	FREE_VECTOR(x);
}

double generate_normal(double mu, double sigma) 
{
	double U, V;
	U = ((float)rand())/RAND_MAX;
	V = ((float)rand())/RAND_MAX;
	return(mu + sigma*sqrt(-2*log(U))*cos(2*M_PI*V));
}

int max_vec(double *vec, int n)
// Returns the index with the largest value out of vec
{
	int i, m=0;
	double max=vec[0];
	for(i=1; i<n; ++i){
		if(vec[i]>max){
			m = i;
			max = vec[i];
		}
	}
	return m;
}

int min(int a, int b)
{
	if(a<=b)	return a;
	else		return b;
}

int max(int a, int b)
{
	if(a>=b)	return a;
	else		return b;
}

double dist_max(double *a, double *b, int n)
{
	int i;
	double max;
	double *diff;		MAKE_VECTOR(diff, n);
	for(i=0; i<n; ++i)	diff[i]=fabs(a[i]-b[i]);
	max = diff[0];
	for(i=1; i<n; ++i){
		if( diff[i] >max)
			max = diff[i];
	}
	return max;
}		

void par_update(int p, int q, double *alpha, double *beta, double sig2, double *par)
{
	int i;
	for(i=0; i<q; ++i)	par[i] = beta[i];
	par[q] = sig2;
	for(i=0; i<p; ++i)	par[q+1+i] = alpha[i];
}

void zero_mat(int nrow, int ncol, double **mat)
{
	int i, j;
	for(i=0; i<nrow; ++i){
		for(j=0; j<ncol; ++j){
			mat[i][j] = 0.;
		}
	}
}

void OLS(double **X, double *y, int n, int q, double *beta)
{
	double **XpX; 	MAKE_MATRIX(XpX, n, q);
	double *Xpy;	MAKE_VECTOR(Xpy, q);
	X_transpose_X(X, n, q, XpX);
	my_inv(q, XpX);
	X_transpose_y(X, y, n, q, Xpy);
	matxvec(XpX, q, q, Xpy, q, beta);
	FREE_MATRIX(XpX);
	FREE_VECTOR(Xpy);
}

double samp_var(double **X, double *r, int n, int q, double *beta)
{
	int t;
	double sum=0.;
	double *mu;		MAKE_VECTOR(mu, n);
	matxvec(X, n, q, beta, q, mu);
	for(t=0; t<n; ++t)		sum += (r[t]-mu[t]) * (r[t]-mu[t]);
	return sum / n;
	FREE_VECTOR(mu);
}

void X_transpose_X(double **X, int nrow, int ncol, double **XpX)
{
	int i, j, k;
	double temp=0;
	for(i=0; i<ncol; ++i){
		for(j=0; j<ncol; ++j){
			for(k=0; k<nrow; ++k){
				temp += X[k][i]* X[k][j];
			}
			XpX[i][j] = XpX[j][i] = temp;
			temp=0;
		}
	}
}

void X_transpose_y(double **X, double *y, int nrow, int ncol, double *Xpy)
{
	int i, k;
	double temp=0;
	for(i=0; i<ncol; ++i){
		for(k=0; k<nrow; ++k){
			temp+=X[k][i]*y[k];
		}
		Xpy[i]=temp;
		temp=0;
	}
}

/* Multiplies matrices a and b and puts the result in c which should be
 * pre-allocated.   exit() will be called if a and b are incompatible
*/
int multiply(double **a, int arows, int acols,
	      double **b, int brows, int bcols, double **c)
{
  int i, j, k;
  
  assert(acols==brows);
  
  for (i=0; i<arows; ++i)
    for (j=0; j<bcols; ++j) {
      c[i][j] = 0;
      for (k=0; k<acols; ++k)
	c[i][j] += a[i][k] * b[k][j];
    }
  return 0;
}

/* Multiplies matrix a and vector x and puts the result in y which should be
 * pre-allocated.   exit() will be called if a and x are incompatible
*/
int matxvec(double **a, int arows, int acols,
		double *x, int xrows, double *y)
{
  int i, k;
  
  assert(acols==xrows);
  
  for (i=0; i<arows; ++i){
    y[i] = 0;
    for (k=0; k<acols; ++k){
      y[i] += a[i][k] * x[k];
    }
  }
  return 0;
}

void transpose(int nrow, int ncol, double **M, double **tM)
{
	int i, j;
	for(i=0; i<nrow; ++i)
		for(j=0; j<ncol; ++j)
			tM[j][i] = M[i][j];
}

/* Prints matrix with a specified format */
void print_dmatrix(double **a, int rows, int cols, const char *format)
{
	int i, j;
	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++)
			printf(format, a[i][j]);
		putchar('\n');
	}
	printf("\n");
}

/* Prints vector with a spefified format */
void print_dvector(double *a, int rows, const char *format)
{
  int i;
  for (i=0; i<rows; i++) {
    printf(format, a[i]);
  }
  printf("\n");
}

double quadratic(double **A,double *x,int p) 
{                  /*calculates x'Ax where A is a p x p-dimensional matrix */
  int i,j;
  double qform=0.0;
  for(i=0;i<p;i++) {
    for(j=0;j<p;j++) {
      qform+=x[i]*x[j]*A[i][j];
    }
  }
  return(qform);
}

/*copy matrix A to matrix B*/
int cpy(double **a,int nrows,int ncols,double **b)
{
  int i,j;
  for(i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      b[i][j]=a[i][j];
    }
  }
  return 0;
}

/* Computes computes the bilinear form t(a)%*%M%*%b */
double bilinear(double *a, double **M, double *b, int n)
{
	int i;
	double ans = 0.0;
	double *Mb;		MAKE_VECTOR(Mb, n);
	
	matxvec(M, n, n, b, n, Mb);
	for(i=0; i<n; ++i)	ans += a[i] * Mb[i];
	return ans;
	FREE_VECTOR(Mb);
}


void my_inv(int size, double **mat)
{
	double det, a, d;
	if(size==1)
		mat[0][0] = 1 / mat[0][0];
	else if(size==2){
		det = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
		mat[0][1] = -1. * mat[0][1] / det;
		mat[1][0] = -1.* mat[1][0] / det;
		a = mat[0][0];	d = mat[1][1];
		mat[0][0] = d / det;
		mat[1][1] = a / det;
	}
	else
		matinv(size, mat, &det);
}

double my_det(int size, double **mat)
{
	double det;
	if(size==1)
		return mat[0][0];
	if(size==2)
		return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
	else{
		matinv(size, mat, &det);
		return det;
	}
}


/* Matrix inverse. Use posymatinv for inversion of a positive-definite 
   matrix. 
  
   Author: Ranjan Maitra <maitra@iastate.edu>
   Date:   03/07/2005
   Uses:   LAPACK 
*/  

#include <string.h>

void dgetrf_(int *Mp, int *Np, double *A, int *LDA, int *PIVOT, int *INFOp);
void dgetri_(int *Np, double *A, int *LDA, int *PIVOT, double *WORK, 
	     int *LWORK, int *INFOp);

void  dpotrf_(char *UPLOp,int *Np, double *A, int *LDAp, int *INFOp);
void  dpotri_(char *UPLOp,int *Np, double *A, int *LDAp, int *INFOp);

void dpptrf_(char *UPLOp,int *Np,double *A,int *INFOp);
void dpptri_(char *UPLOp,int *Np,double *A,int *INFOp);

int matinv(int sizeA,double **A,double (*determinant))
{
  int i, j , *pivot,N=sizeA*sizeA,size=sizeA;
  double *AT,*work;	/* AT=transpose vectorized matrix (to accomodate
			   Fortran) 
			   work=workspace vector */
  int INFO,ipiv=1;

  MAKE_VECTOR(AT,size*size);
  MAKE_VECTOR(work,size*size);
  MAKE_VECTOR(pivot,size);

  for (i=0; i<size; i++)		/* to call a Fortran routine from C */
    {				/* have to transform the matrix */
      for(j=0; j<size; j++) AT[j+size*i]=A[j][i];		
    }						
  
  dgetrf_(&size,&size,AT,&size,pivot,&INFO);
  /* LAPACK routine DGETRF computes an LU factorization of a general 
     m x n matrix A using partial pivoting with row interchanges. The
     factorization has the form A = P * L * U where P is a permutation
     matrix, L is lower triangular with unit diagonal elements (lower 
     trapezoidal if m > n), and U is upper triangular (upper trapezoidal
     if m < n). Note that because of the permutation, the determinant 
     needs to be multiplied by -1 for every interchange that has occurred.
     Parameters in the order as they appear in the function call: 
     number of rows of the matrix A, number of columns of the 
     matrix A, the matrix A, the leading dimension of A, the 
     array that records pivoting, and the flag for the
     result. On exit, A contains the factors of L and U (with the
     diagonals of L not stored).*/	  
  if (INFO==0) {
    for(i=0;i<size;i++) {
      if (i!=(pivot[i]-1)) ipiv*=-1; /* PIVOT assumes indices are from 1 
					through N*/
    }
    (*determinant)=(double)ipiv;
    for (i=0;i<size;i++) {
      (*determinant)*=AT[i+i*size];
    }
    dgetri_(&size,AT,&size,pivot,work,&N,&INFO); 
    /* LAPACK routine DGETRI computes the inverse of a matrix A 
       using the output of DGETRF. This method inverts U and then
       computes A^(-1) by solving A^(-1)L = U^(-1) for A^(-1).
       parameters in the order as they appear in the function call:
       order of the matrix A, the matrix A, the leading dimension of
       A, the array that records pivoting, workspace, the 
       dimension of the workspace array, and the flag for the 
       result. On exit, A contains the inverted matrix. */
    if (INFO!=0) {
      printf("Problem in matinv: dgetri error %d\n",INFO);
    }
  }
  else {
    printf("Problem in matinv: dgetrf error %d\n",INFO);
  }
  for (i=0; i<size; i++)		/* to call a Fortran routine from C */
    {				/* have to transform the matrix */
      for(j=0; j<size; j++) {
	A[j][i]=AT[j+size*i];		
      }
    }
  FREE_VECTOR(AT);
  FREE_VECTOR(pivot);
  FREE_VECTOR(work);
  return 0;
}

	

