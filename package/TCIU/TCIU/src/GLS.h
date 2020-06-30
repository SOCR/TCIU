#ifndef GLS_H
#define GLS_H

void GLS(int n, int q, int p, double **X, double **Rninv, double *y, 
	double *beta);

void comp_XRX(int n, int q, int p, double **X, double **Rninv, double **XRX);

void comp_XRy(int n, int q, int p, double **X, double **Rninv, double *y, 
	double *XRy);

double comp_aRb(int n, int p, double **Rninv, double *a, double *b);

void comp_Rpinv(int p, double *alpha, double **Rpinv);

void comp_Rninv(int n, int p, double *alpha, double **Rninv);

#endif
