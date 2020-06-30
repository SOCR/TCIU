#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include "array.h"
#include "mat_vec.h"
#include "RC_interface.h"
#include "complex_smoother.h"

void Rwrapper_complex_running_line(int *N, int *k, double *r, double *phi, 
	double *tol, double *rho_vec, double *theta_vec, int *n_iter, 
	int *max_iter, int *ma_line, int *interp)
{
	int t=0, i, start;
	int m = *N / *interp;
	if(*ma_line==1){ //running line smoother
		for(i=0; i<m; ++i){
			smooth_1t(t, *N, 2, *k, r, phi, *tol, &rho_vec[t], 
				&theta_vec[t], &n_iter[t], *max_iter);
			t = t + *interp;
		}
		start = (m-1)* (*interp);
		for(t=start+1; t<*N; ++t)
			smooth_1t(t, *N, 2, *k, r, phi, *tol, &rho_vec[t], 
				&theta_vec[t], &n_iter[t], *max_iter);
	}
	if(*ma_line==0){ //moving average smoother
		for(t=0; t<*N; ++t){
			ma_1t(t, *N, *k, r, phi, &rho_vec[t], &theta_vec[t]);
	}}
}

void smooth_1t(int t, int N, int q, int k, double *r, double *phi, double tol, 
	double *rho, double *theta, int *n_iter, int max_iter)
{
	int i;
	double gamma, gam0;
	double *beta;		MAKE_VECTOR(beta, q);
	int nhood_min = max(t - k, 0);
	int nhood_max = min(t + k, N-1);
	int n_t = nhood_max - nhood_min + 1;
	double t_c = (double) t - (nhood_min + nhood_max) / 2.0;
	double *r_t;		MAKE_VECTOR(r_t, n_t);
	double *phi_t;		MAKE_VECTOR(phi_t, n_t);
	double **X;			MAKE_MATRIX(X, n_t, q);
	double *u;			MAKE_VECTOR(u, n_t);
	for(i=0; i<n_t; ++i){
		r_t[i] = r[nhood_min + i];
		phi_t[i] = phi[nhood_min + i];
		X[i][0] = 1.0;
		X[i][1] = u[i] = (double) i - (n_t - 1) / 2.0;
	}
	est_par(n_t, q, u, X, &gamma, &gam0, beta, r_t, phi_t, tol, n_iter, 
		max_iter);
	*rho = beta[0] + t_c * beta[1];
	*theta = gam0 + 2.0*atan(t_c*gamma);
	FREE_VECTOR(beta);
	FREE_VECTOR(r_t);			FREE_VECTOR(phi_t);
	FREE_MATRIX(X);				FREE_VECTOR(u);
}

void est_par(int n, int q, double *u, double **X, double *gamma, double *gam0, 
	double *beta, double *r, double *phi, double tol, int *n_iter, int max_iter)
{
	int iter=0, conv=0;
	double LL_old, LL_new;
	starting_values(n, q, r, phi, gam0, gamma, beta);
	while(conv == 0 && iter <= max_iter){
		if(iter==0){
			compute_LL(n, q, u, X, *gamma, *gam0, beta, r, phi, &LL_old);
				// in phase_model.c
		}
		else
			LL_old = LL_new;
		
		++iter;
		one_iter(n, q, u, X, gamma, gam0, beta, r, phi);
		compute_LL(n, q, u, X, *gamma, *gam0, beta, r, phi, &LL_new);
		if(LL_new - LL_old < tol) conv = 1;
	}
	if(iter > max_iter) Rprintf("Over max_iter - complex\n");
	*n_iter = iter;
}

void starting_values(int n, int q, double *r, double *phi, double *gam0, 
	double *gamma, double *beta)
{
	int t, i;
	double s=0, c=0, temp=0;
	for(t=0; t<n; ++t){
		s += r[t] * sin(phi[t]);
		c += r[t] * cos(phi[t]);
	}
	*gam0 = atan2(s, c);
	*gamma = 0.0;
	for(t=0; t<n; ++t) temp += r[t] * cos(phi[t] - *gam0);
	beta[0] = temp / n;
	for(i=1; i<q; ++i) beta[i] = 0.0;
}

void compute_LL(int n, int q, double *u, double **X, double gamma, double gam0, 
	double *beta, double *r, double *phi, double *LL)
{
	int t;
	double *rho;		MAKE_VECTOR(rho, n);
	double temp = 0.0, sig2;
	matxvec(X, n, q, beta, q, rho);
	for(t=0; t<n; ++t)
		temp += r[t]*r[t] + rho[t]*rho[t] - 
				2.0*rho[t]*r[t] * cos(phi[t] - gam0 - 2.0*atan(u[t]*gamma));
	sig2 = temp / (2.0 * n);
	*LL = -1.0 * n * log(sig2) - n;
	FREE_VECTOR(rho);
}

void one_iter(int n, int q, double *u, double **X, double *gamma, double *gam0, 
	double *beta, double *r, double *phi)
{
	int t;
	double *g;			MAKE_VECTOR(g, n);
	double *g1;			MAKE_VECTOR(g1, n);
	double *g2;			MAKE_VECTOR(g2, n);
	double *phi_res;	MAKE_VECTOR(phi_res, n);
	double *r_star;		MAKE_VECTOR(r_star, n);
	double *rho;		MAKE_VECTOR(rho, n);
	double s=0, c=0, temp, num=0, den=0;
	double gamma_old = *gamma;
	matxvec(X, n, q, beta, n, rho);
	for(t=0; t<n; ++t){
		g[t] = 2.0*atan(u[t] * (*gamma));
		s += r[t] * rho[t] * sin(phi[t] - g[t]);
		c += r[t] * rho[t] * cos(phi[t] - g[t]);
	}
	*gam0 = atan2(s, c);
	for(t=0; t<n; ++t){ 
		phi_res[t] = phi[t] - *gam0 - g[t];
		r_star[t] = r[t] * cos(phi_res[t]);
	}
	OLS(X, r_star, n, q, beta);
	matxvec(X, n, q, beta, n, rho);
	for(t=0; t<n; ++t){ 
		temp = u[t] * (*gamma);
		g1[t] = 2 / (1 + temp*temp);
		g2[t] = 2 * temp / pow(1 + temp*temp, 2);
		num += u[t]*r[t]*rho[t]*g1[t]*sin(phi_res[t]);
		den += u[t]*u[t] * r[t] * rho[t] * (g1[t]*g1[t] * cos(phi_res[t]) - 
											 g2[t] * sin(phi_res[t]));
	}
	*gamma = gamma_old + num / den;
	FREE_VECTOR(g);			FREE_VECTOR(g1);
	FREE_VECTOR(g2);		FREE_VECTOR(phi_res);
	FREE_VECTOR(r_star);	FREE_VECTOR(rho);
}

void Rwrapper_mag_only_run_line(int *N, int *k, double *r, double *rho_vec, 
	int *interp)
{
	int t=0, i, start;
	int m = *N / *interp;
	//Rprintf("m=%d\n", m);
	for(i=0; i<m; ++i){
		mag_smooth_1t(t, *N, *k, r, &rho_vec[t]);
		t = t + *interp;
	}
	start = (m-1)* (*interp);
	for(i=start+1; i<*N; ++i)
		mag_smooth_1t(i, *N, *k, r, &rho_vec[i]);
}

void mag_smooth_1t(int t, int N, int k, double *r, double *rho)
{
	int i;
	double *beta;		MAKE_VECTOR(beta, 2);
	int nhood_min = max(t - k, 0);
	int nhood_max = min(t + k, N-1);
	int n_t = nhood_max - nhood_min + 1;
	double t_c = (double) t - (nhood_min + nhood_max) / 2.0;
	double *r_t;		MAKE_VECTOR(r_t, n_t);
	double **X;			MAKE_MATRIX(X, n_t, 2);
	for(i=0; i<n_t; ++i){
		r_t[i] = r[nhood_min + i];
		X[i][0] = 1.0;
		X[i][1] = (double) i - (n_t - 1) / 2.0;
	}
	OLS(X, r_t, n_t, 2, beta);
	*rho = beta[0] + t_c * beta[1];
	FREE_VECTOR(beta);
	FREE_VECTOR(r_t);
	FREE_MATRIX(X);
}

void ma_1t(int t, int N, int k, double *r, double *phi,  
	double *rho, double *theta)
{
	int i;
	int nhood_min = max(t - k, 0);
	int nhood_max = min(t + k, N-1);
	int n_t = nhood_max - nhood_min + 1;
	double *r_t;		MAKE_VECTOR(r_t, n_t);
	double *phi_t;		MAKE_VECTOR(phi_t, n_t);
	for(i=0; i<n_t; ++i){
		r_t[i] = r[nhood_min + i];
		phi_t[i] = phi[nhood_min + i];
	}
	ma(n_t, r_t, phi_t, rho, theta);
	FREE_VECTOR(r_t); FREE_VECTOR(phi_t);
}

void ma(int n, double *r, double *phi, double *rho, double *theta)
{
	int t;
	double s=0, c=0, temp=0;
	for(t=0; t<n; ++t){
		s += r[t] * sin(phi[t]);
		c += r[t] * cos(phi[t]);
	}
	*theta = atan2(s, c);
	for(t=0; t<n; ++t) temp += r[t] * cos(phi[t] - *theta);
	*rho = temp / n;
}	




/*void res_beta(int n, double **X, double *beta, double *r, double *phi, 
	double *u, double gam0, double gamma)
{
	int i;
	double **C;				MAKE_MATRIX(C, 1, 2);
	double **Psi;			MAKE_MATRIX(Psi, 2, 2);
	double **Rninv;			MAKE_MATRIX(Rninv, 1, 1); //not used
	double *betahat;		MAKE_VECTOR(betahat, 2);
	double *beta0;			MAKE_VECTOR(beta0, 2);
	double *betan;			MAKE_VECTOR(betan, 2);
	double ss0, ssn;
	double mu0=0, mun=0;
	for(i=0; i<2; ++i){
		betahat[i] = beta[i];
		mu0 += X[0][i] * betahat[i];
		mun += X[n-1][i] * betahat[i];
	}
	if(mu0 < 0.0){
		C[0][0] = X[0][0];	C[0][1] = X[0][1];
		comp_Psi(n, 2, 0, 1, X, Rninv, 1, C, Psi);
		matxvec(Psi, 2, 2, betahat, 2, beta0);
		for(i=0; i<2; ++i) beta[i] = beta0[i];
	}
	if(mun < 0.0){
		C[0][0] = X[n-1][0];	C[0][1] = X[n-1][1];
		comp_Psi(n, 2, 0, 1, X, Rninv, 1, C, Psi);
		matxvec(Psi, 2, 2, betahat, 2, betan);
		for(i=0; i<2; ++i) beta[i] = betan[i];
	}
	if(mu0 < 0.0 && mun < 0.0){
		ss0 = compute_ss(n, X, beta0, r, phi, u, gam0, gamma);
		ssn = compute_ss(n, X, betan, r, phi, u, gam0, gamma);
		if(ss0 < ssn) 
			for(i=0; i<2; ++i) beta[i] = beta0[i];
		else 
			for(i=0; i<2; ++i) beta[i] = betan[i];
	}
	//Rprintf("\n");
	FREE_MATRIX(C);			FREE_MATRIX(Psi);
	FREE_MATRIX(Rninv);		FREE_VECTOR(betahat);
	FREE_VECTOR(beta0);		FREE_VECTOR(betan);
}

double compute_ss(int n, double **X, double *beta, double *r, double *phi, 
	double *u, double gam0, double gamma)
{
	int t;
	double *rho;		MAKE_VECTOR(rho, n);
	double *phi_res;	MAKE_VECTOR(phi_res, n);
	double ss = 0;
	matxvec(X, n, 2, beta, 2, rho);
	for(t=0; t<n; ++t){
		phi_res[t] = phi[t] - gam0 - 2.0 * atan(u[t]*gamma);
		ss += r[t]*r[t] + rho[t]*rho[t] - 2.0 * rho[t]*r[t]*cos(phi_res[t]);
	}
	FREE_VECTOR(rho);
	FREE_VECTOR(phi_res);
	return ss;
}*/