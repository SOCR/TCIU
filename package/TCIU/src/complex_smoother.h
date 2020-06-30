#ifndef COMPLEX_SMOOTHER_H
#define COMPLEX_SMOOTHER_H

void Rwrapper_complex_running_line(int *N, int *k, double *r, double *phi, 
	double *tol, double *rho_vec, double *theta_vec, int *n_iter, 
	int *max_iter, int *ma_line, int *interp);

void Rwrapper_mag_only_run_line(int *N, int *k, double *r, double *rho_vec, 
	int *interp);

void mag_smooth_1t(int t, int N, int k, double *r, double *rho);

void smooth_1t(int t, int N, int q, int k, double *r, double *phi, double tol, 
	double *rho, double *theta, int *n_iter, int max_iter);

void est_par(int n, int q, double *u, double **X, double *gamma, double *gam0, 
	double *beta, double *r, double *phi, double tol, int *n_iter, 
	int max_iter);

void starting_values(int n, int q, double *r, double *phi, double *gam0, 
	double *gamma, double *beta);

void compute_LL(int n, int q, double *u, double **X, double gamma, double gam0, 
	double *beta, double *r, double *phi, double *LL);

void ma_1t(int t, int N, int k, double *r, double *phi,  
	double *rho, double *theta);

void ma(int n, double *r, double *phi, double *rho, double *theta);

void one_iter(int n, int q, double *u, double **X, double *gamma, double *gam0, 
	double *beta, double *r, double *phi);

void res_beta(int n, double **X, double *beta, double *r, double *phi, 
	double *u, double gam0, double gamma);

double compute_ss(int n, double **X, double *beta, double *r, double *phi, 
	double *u, double gam0, double gamma);

#endif
