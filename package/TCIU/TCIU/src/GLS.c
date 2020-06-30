#include <stdio.h>		/* for fprintf() */
#include <stdlib.h>		/* for malloc() */
#include <math.h>		/* for sqrt() */
#include <assert.h>
#include "array.h"
#include "mat_vec.h"
#include "GLS.h"

void GLS(int n, int q, int p, double **X, double **Rninv, double *y, 
	double *beta)
{
	double **XRX;		MAKE_MATRIX(XRX, q, q);
	double *XRy;		MAKE_VECTOR(XRy, q);
	comp_XRX(n, q, p, X, Rninv, XRX);
	comp_XRy(n, q, p, X, Rninv, y, XRy);
	my_inv(q, XRX);
	matxvec(XRX, q, q, XRy, q, beta);
	FREE_MATRIX(XRX);
	FREE_VECTOR(XRy);
}

void comp_XRX(int n, int q, int p, double **X, double **Rninv, double **XRX)
{
	int i, j, t, tp;
	double temp = 0.;
	for(i=0; i<q; ++i){
		for(j=i; j<q; ++j){
			for(t=0; t<n; ++t){
				for(tp=max(t-p,0); tp<=min(t+p,n-1); ++tp)
					temp += X[t][i]*Rninv[t][tp]*X[tp][j];
			}		
			XRX[i][j] = XRX[j][i] = temp;
			temp = 0.;
		}
	}
}

void comp_XRy(int n, int q, int p, double **X, double **Rninv, double *y, 
	double *XRy)
{
	int i, t, tp;
	double temp=0.;
	for(i=0; i<q; ++i){
		for(t=0; t<n; ++t){
			for(tp=max(t-p,0); tp<=min(t+p,n-1); ++tp)
				temp += X[t][i]*Rninv[t][tp]*y[tp];
		}
		XRy[i] = temp;
		temp = 0.;
	}
}

double comp_aRb(int n, int p, double **Rninv, double *a, double *b)
//computes bilinear form t(a) %*% Rninv %*% b
{
	int t, tp;
	double ans = 0.0;
	for(t=0; t<n; ++t)
		for(tp=max(t-p,0); tp<=min(t+p,n-1); ++tp)
			ans += a[t] * Rninv[t][tp] * b[tp];
	return ans;
}

void comp_Rninv(int n, int p, double *alpha, double **Rninv)
{
	int i, j, k;
	double temp = 0.0;
	double *atil;		MAKE_VECTOR(atil, p+1);
	atil[0] = 1.0;	
	for(i=1; i<=p; ++i)	atil[i] = -1.0*alpha[i-1];
	//Fills corners of Rninv
	for(k=0; k<=p-1; ++k){ //abs(row-col)
		for(i=0; i<=p-1-k; ++i){ //row, col
			for(j=0; j<=i; ++j)	temp+= atil[j]*atil[j+k];
			Rninv[i][i+k] = Rninv[i+k][i] = Rninv[n-i-k-1][n-i-1] = 
				Rninv[n-i-1][n-i-k-1] = temp;
			temp = 0.0;
		}
	}
	//Fills rest of Rn
	for(k=0; k<=p; ++k){
		for(j=0; j<=p-k; ++j)	temp+= atil[j]*atil[j+k];
		for(i = p; i<n-p; ++i){
			Rninv[i][i+k] = Rninv[i+k][i] = temp;
		}
		if(k>0)
			for(i=0; i<=k-1; ++i)
				Rninv[p+i][p-k+i] = 
					Rninv[p-k+i][p+i] = temp;
		temp = 0.0;
	}
	FREE_VECTOR(atil);
}

void comp_Rpinv(int p, double *alpha, double **Rpinv)
{
	int i, j, k;
	double temp = 0.;
	double **Rn;		MAKE_MATRIX(Rn, p, p);
	double *atil;			MAKE_VECTOR(atil, p+1);
	atil[0] = 1.0;
	for(i=0; i<p; ++i)	atil[i+1] = -1.0 * alpha[i];
	//Fills corner of Rn
	for(k=0; k<=p-1; ++k){
		for(i=0; i<=p-1-k; ++i){
			for(j=0; j<=i; ++j)	temp += atil[j]*atil[j+k];
			if(k==0)
				Rn[i][i] = temp;
			else
				Rn[i][i+k] = Rn[i+k][i] = temp;
			temp = 0.;
		}
	}
	//Creates Rpinv
	for(i=0; i<=p-1; ++i){
		for(k=0; k<=p-i-1; ++k){
			for(j=0; j<=i; ++j)	temp += atil[p-j]*atil[p-j-k];
			Rpinv[i][i+k] = Rpinv[i+k][i] =Rn[i][i+k] - temp;	
			 temp = 0.;
		}
	}
	FREE_MATRIX(Rn);
	FREE_VECTOR(atil);
}