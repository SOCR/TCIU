#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "array.h"
#include "mat_vec.h"
#include "RC_interface.h"
#include "GLS.h"
#include "complex_Sig_sig2I.h"

void Rwrapper_complex_activation(int *n, int *q, int *p, int *m, 
	double *C_vec, double *X_vec, double *yR, double *yI, 
	int *max_iter, double *LL_eps, double *par_vec, 
	double *lrt_stat)
{
	int len = *q + 2 + *p, i;
	double **X;			MAKE_MATRIX(X, *n, *q);
	double **C;			MAKE_MATRIX(C, *m, *q);
	double **par;			MAKE_MATRIX(par, 2, len);
	double *par_unres;		MAKE_VECTOR(par_unres, len);
	double *par_res;		MAKE_VECTOR(par_res, len);
	copy_1d_to_2d(X_vec, X, *n, *q);
	copy_1d_to_2d(C_vec, C, *m, *q);
	complex_activation(*n, *q, *p, *m, C, X, yR, yI, *max_iter, 
		*LL_eps, par_unres, par_res, lrt_stat);
	for(i=0; i<len; ++i){
		par[0][i] = par_unres[i];
		par[1][i] = par_res[i];
	}
	copy_2d_to_1d(par, 2, len, par_vec);
	FREE_MATRIX(par);
	FREE_VECTOR(par_unres);
	FREE_VECTOR(par_res);
	FREE_MATRIX(X);
	FREE_MATRIX(C);
}

void complex_activation(int n, int q, int p, int m, double **C, 
	double **X, double *yR, double *yI, int max_iter, 
	double LL_eps, double *par_unres, double *par_res, double *lrt_stat)
{
	if(p==0)
		complex_indep(n, q, m, C, X, yR, yI, par_unres, par_res, lrt_stat);
	else 
		complex_arp(n, q, p, m, X, C, yR, yI, max_iter, 
			LL_eps, par_unres, par_res, lrt_stat);
}

//Independent functions

void complex_indep(int n, int q, int m, double **C, double **X, 
	double *yR, double *yI, double *par_unres, double *par_res,
	double *lrt_stat)
{
	int i;
	double *beta_hat;		MAKE_VECTOR(beta_hat, q);
	double *beta_til;		MAKE_VECTOR(beta_til, q);
	double sig2_hat, theta_hat, sig2_til, theta_til;
	complex_indep_unres(X, yR, yI, n, q, beta_hat, &sig2_hat, &theta_hat);
	complex_res_indep(n, q, m, C, X, yR, yI, beta_til, &sig2_til, &theta_til);
	//copy unres pars
	for(i=0; i<q; ++i)	par_unres[i] = beta_hat[i];
	par_unres[q] = theta_hat;
	par_unres[q+1] = sig2_hat;
	//copy res pars
	for(i=0; i<q; ++i)	par_res[i] = beta_til[i];
	par_res[q] = theta_til;
	par_res[q+1] = sig2_til;
	*lrt_stat = 2.0 * n * log(sig2_til / sig2_hat);
	FREE_VECTOR(beta_hat);
	FREE_VECTOR(beta_til);
}

void complex_indep_unres(double **X, double *yR, double *yI, int n, int q, 
	double *beta, double *sig2, double *theta)
{
	int t, i;
	double num, den, temp = 0.;
	double *bR;		MAKE_VECTOR(bR, q);
	double *bI;		MAKE_VECTOR(bI, q);
	double **XpX;		MAKE_MATRIX(XpX, q, q); 
	double *mu;		MAKE_VECTOR(mu, n);
	OLS(X, yR, n, q, bR);
	OLS(X, yI, n, q, bI);
	X_transpose_X(X, n, q, XpX);
	num = 2.*bilinear(bR, XpX, bI, q);
	den = quadratic(XpX, bR, q) - quadratic(XpX, bI, q);
	*theta = 0.5*atan2(num, den);
	for(i=0; i<q; ++i)	beta[i] = bR[i]*cos(*theta)+bI[i]*sin(*theta);
	matxvec(X, n, q, beta, q, mu);
	for(t=0; t<n; ++t)	
		temp+=(yR[t]-mu[t]*cos(*theta))*(yR[t]-mu[t]*cos(*theta))+
				(yI[t]-mu[t]*sin(*theta))*(yI[t]-mu[t]*sin(*theta));
	*sig2 = temp/(2.0*n);
	//identifiability issue
	if(beta[0]<0){
		for(i=0; i<q; ++i)	beta[i] = -1.*beta[i];
		if(*theta>0)		
			*theta = *theta - M_PI;
		else
			*theta = *theta +M_PI;
	}
	FREE_VECTOR(bR);
	FREE_VECTOR(bI);
	FREE_MATRIX(XpX);
	FREE_VECTOR(mu);
}

void complex_res_indep(int n, int q, int m, double **C, double **X, 
	double *yR, double *yI, double *beta, double *sig2, double *theta)
{
	int t, i;
	double num, den, temp = 0.;
	double *bR;		MAKE_VECTOR(bR, q);
	double *bI;		MAKE_VECTOR(bI, q);
	double **Psi;		MAKE_MATRIX(Psi, q, q);
	double *Psi_bR;	MAKE_VECTOR(Psi_bR, q);
	double *Psi_bI;	MAKE_VECTOR(Psi_bI, q);
	double **XpX;		MAKE_MATRIX(XpX, q, q); 
	double *mu;		MAKE_VECTOR(mu, n);
	double **I_n;		MAKE_MATRIX(I_n, n, n);
	make_identity_mat(n, I_n);
	OLS(X, yR, n, q, bR);
	OLS(X, yI, n, q, bI);
	comp_Psi(n, q, 0, m, X, I_n, 1, C, Psi);
	matxvec(Psi, q, q, bR, q, Psi_bR);
	matxvec(Psi, q, q, bI, q, Psi_bI);
	X_transpose_X(X, n, q, XpX);
	num = 2.*bilinear(Psi_bR, XpX, bI, q);
	den = bilinear(Psi_bR, XpX, bR, q) - bilinear(Psi_bI, XpX, bI, q);
	*theta = 0.5 * atan2(num, den);
	for(i=0; i<q; ++i)	beta[i] = Psi_bR[i]*cos(*theta) + Psi_bI[i]*sin(*theta);
	matxvec(X, n, q, beta, q, mu);
	for(t=0; t<n; ++t)	
		temp+=(yR[t]-mu[t]*cos(*theta))*(yR[t]-mu[t]*cos(*theta))+
			(yI[t]-mu[t]*sin(*theta))*(yI[t]-mu[t]*sin(*theta));
	*sig2 = temp/(2.0*n);
	//identifiability issue
	if(beta[0]<0){
		for(i=0; i<q; ++i)	beta[i] = -1.0 * beta[i];
		if(*theta>0)	*theta = *theta - M_PI;
		else 		*theta = *theta +M_PI;
	}
	FREE_VECTOR(bR);
	FREE_VECTOR(bI);
	FREE_MATRIX(Psi);
	FREE_VECTOR(Psi_bR);
	FREE_VECTOR(Psi_bI);
	FREE_MATRIX(XpX);
	FREE_VECTOR(mu);
	FREE_MATRIX(I_n);
}

void comp_Psi(int n, int q, int p, int m, double **X, double **Rninv, 
	int Rn_is_In, double **C, double **Psi)
{
	int i, j;
	double **XRX;		MAKE_MATRIX(XRX, q, q);
	double **Iq;		MAKE_MATRIX(Iq, q, q);
	double **tC;		MAKE_MATRIX(tC, q, m);
	double **temp1;	MAKE_MATRIX(temp1, q, m);
	double **CXRX;	MAKE_MATRIX(CXRX, m, q);
	double **CXRXC;	MAKE_MATRIX(CXRXC, m, m);
	double **temp2;	MAKE_MATRIX(temp2, q, m);
	double **temp3;	MAKE_MATRIX(temp3, q, q);
	if(Rn_is_In == 1) X_transpose_X(X, n, q, XRX);
	else comp_XRX(n, q, p, X, Rninv, XRX);
	my_inv(q, XRX);		//now XRX is XRXinv
	make_identity_mat(q, Iq);
	transpose(m, q, C, tC);
	multiply(XRX, q, q, tC, q, m, temp1);
	multiply(C, m, q, XRX, q, q, CXRX);
	multiply(CXRX, m, q, tC, q, m, CXRXC);
	my_inv(m, CXRXC);		//now CXRXC is (CXRXC)^-1
	multiply(temp1, q, m, CXRXC, m, m, temp2);
	multiply(temp2, q, m, C, m, q, temp3);
	for(i=0; i<q; ++i)
		for(j=0; j<q; ++j)
			Psi[i][j] = Iq[i][j] - temp3[i][j];
	FREE_MATRIX(XRX);
	FREE_MATRIX(Iq);
	FREE_MATRIX(tC);
	FREE_MATRIX(temp1);
	FREE_MATRIX(CXRX);
	FREE_MATRIX(CXRXC);
	FREE_MATRIX(temp2);
	FREE_MATRIX(temp3);
}

// AR(p) functions

void complex_arp(int n, int q, int p, int m, double **X, double **C, 
	double *yR, double *yI, int max_iter, double LL_eps, 
	double *par_unres, double *par_res, double *lrt_stat)
{
	complex_unres_est(n, q, p, X, yR, yI, max_iter,
		par_unres, LL_eps);
	complex_res_est(n, q, p, m, X, C, yR, yI, max_iter,
		par_res, LL_eps);
	compute_lrt_stat(n, p, q, 0, par_unres, par_res, lrt_stat);
}

void complex_unres_est(int n, int q, int p, double **X, 
	double *yR, double *yI, int max_iter,
	double *par, double LL_eps)
{
	int iter=0, conv=0;
	double **D;		MAKE_MATRIX(D, p+1, p+1);
	double **Rninv;	MAKE_MATRIX(Rninv, n, n);
	double *beta;		MAKE_VECTOR(beta, q);
	double *alpha;	MAKE_VECTOR(alpha, p);
	double theta, sig2, LL_old, LL_new;
	make_identity_mat(n, Rninv);		//Rninv = I_n
	comp_theta_beta_arp(n, q, p, X, yR, yI, Rninv, beta, &theta);
	complex_compD(n, p, q, X, beta, theta, yR, yI, D);
	complex_comp_alpha(n, p, D, alpha);
	comp_Rninv(n, p, alpha, Rninv); // in GLS.c
	complex_comp_sig2(n, p, q, yR, yI, X, beta, theta, Rninv, &sig2);
	LL_old = compute_logL(p, n, alpha, sig2, 0);
	while(conv==0 && iter<max_iter){
		++iter;
		comp_theta_beta_arp(n, q, p, X, yR, yI, Rninv, beta, &theta);
		complex_compD(n, p, q, X, beta, theta, yR, yI, D);
		complex_comp_alpha(n, p, D, alpha);
		comp_Rninv(n, p, alpha, Rninv);
		complex_comp_sig2(n, p, q, yR, yI, X, beta, theta, Rninv, &sig2);
		LL_new = compute_logL(p, n, alpha, sig2, 0);
		if(LL_new - LL_old < LL_eps) conv =1;
		else LL_old = LL_new;
	}
	store_est(q, p, beta, theta, sig2, alpha, par);
	FREE_MATRIX(D);
	FREE_MATRIX(Rninv);
	FREE_VECTOR(beta);
	FREE_VECTOR(alpha);	
}

void comp_theta_beta_arp(int n, int q, int p, double **X, 
	double *yr, double *yi, double **Rninv,
	double *beta, double *theta)
{
	int i;
	double num, den;
	double **XRX;		MAKE_MATRIX(XRX, q, q);
	double *beta_r;		MAKE_VECTOR(beta_r, q);
	double *beta_i;		MAKE_VECTOR(beta_i, q);
	GLS(n, q, p, X, Rninv, yr, beta_r); // in GLS.c
	GLS(n, q, p, X, Rninv, yi, beta_i);
	comp_XRX(n, q, p, X, Rninv, XRX); // in GLS.c
	num = 2.*bilinear(beta_r, XRX, beta_i, q);
	den = quadratic(XRX, beta_r, q) - quadratic(XRX, beta_i, q);
	*theta= 0.5*atan2(num, den);
	for(i=0; i<q; ++i)	beta[i] = beta_r[i]*cos(*theta)+beta_i[i]*sin(*theta);
	//identifiability issue
	if(beta[0]<0){
		for(i=0; i<q; ++i)	beta[i] = -1.*beta[i];
		if(*theta>0)		
			*theta = *theta - M_PI;
		else
			*theta = *theta +M_PI;
	}
	FREE_MATRIX(XRX);
	FREE_VECTOR(beta_r);
	FREE_VECTOR(beta_i);
}

void complex_compD(int n, int p, int q, double **X, double *beta, double theta, double *yr,
	double *yi, double **D)
{
	int t, i, k;
	double temp=0.0, cosine=cos(theta), sine=sin(theta);
	double *mu;		MAKE_VECTOR(mu, n);
	matxvec(X, n, q, beta, q, mu);
	for(i=0; i<=p; ++i){
		for(k=0; k<=p-i; ++k){
			for(t=i; t<n-i-k; ++t){
				temp+= (yr[t]-mu[t]*cosine)*(yr[t+k]-mu[t+k]*cosine) +
						(yi[t]-mu[t]*sine)*(yi[t+k]-mu[t+k]*sine);
			}
			D[i][i+k] = D[i+k][i] = temp;
			temp=0.0;
		}
	}
	FREE_VECTOR(mu);
}

void complex_comp_alpha(int n, int p, double **D, double *alpha)
{
	int i, j;
	double *gamma;	MAKE_VECTOR(gamma, p+1);
	double **alp_mat;	MAKE_MATRIX(alp_mat, p, p);
	double *d;			MAKE_VECTOR(d, p);
	//estimate gamma from D
	for(i=0; i<=p; ++i)	gamma[i] = D[i][0] / (2.*n);
	//solve for alpha
	for(i=1; i<=p; ++i)	d[i-1] = D[i][0];
	for(i=1; i<=p; ++i)
		for(j=1; j<=p; ++j)		
			alp_mat[i-1][j-1] = D[i][j] + 2.*j*gamma[abs(j-i)];
	my_inv(p, alp_mat);
	matxvec(alp_mat, p, p, d, p, alpha);
	FREE_VECTOR(gamma);
	FREE_MATRIX(alp_mat);
	FREE_VECTOR(d);
}

void complex_comp_sig2(int n, int p, int q, double *yR, double *yI, double **X, double *beta, 
	double theta, double **Rninv, double *sig2)
{
	int t, tp;
	double *mu;		MAKE_VECTOR(mu, n);
	double *eR;		MAKE_VECTOR(eR, n);
	double *eI;		MAKE_VECTOR(eI, n);
	double temp=0., cosine=cos(theta), sine=sin(theta);
	matxvec(X, n, q, beta, q, mu);
	for(t=0; t<n; ++t){
		eR[t] = yR[t] - mu[t] * cosine;
		eI[t] = yI[t] - mu[t] * sine;
	}
	for(t=0; t<n; ++t){
		for(tp=max(t-p,0); tp<=min(t+p,n-1); ++tp){
			temp += eR[t]*Rninv[t][tp]*eR[tp] + eI[t]*Rninv[t][tp]*eI[tp];
		}
	}
	*sig2 = temp / (2.0*n);
	FREE_VECTOR(mu);
	FREE_VECTOR(eR);
	FREE_VECTOR(eI);
}

double compute_logL(int p, int n, double *alpha, double sig2, 
	int complex_mag)
{
	double **Rpinv;		MAKE_MATRIX(Rpinv, p, p);
	double det;
	comp_Rpinv(p, alpha, Rpinv); //in GLS.c
	det = fabs(my_det(p, Rpinv));
	FREE_MATRIX(Rpinv);
	if(complex_mag==0) //complex
		return -1.0* n * log(sig2) + log(det) - n;
	if(complex_mag==1) //mag
		return -0.5* n * log(sig2) + 0.5 * log(det) - n/2.0;

	return 0.1; // unreachable	
}

void complex_res_est(int n, int q, int p, int m, double **X, 
	double **C, double *yR, double *yI, int max_iter,
	double *par, double LL_eps)
{
	int iter=0, conv=0;
	double **D;		MAKE_MATRIX(D, p+1, p+1);
	double **Rninv;	MAKE_MATRIX(Rninv, n, n);
	double *beta;		MAKE_VECTOR(beta, q);
	double *alpha;	MAKE_VECTOR(alpha, p);
	double theta, sig2, LL_old, LL_new;
	make_identity_mat(n, Rninv);		//Rninv = I_n
	complex_res_beta_theta(n, q, p, m, C, X, yR, yI, beta, 
		&theta, 1, Rninv);
	complex_compD(n, p, q, X, beta, theta, yR, yI, D);
	complex_comp_alpha(n, p, D, alpha);
	comp_Rninv(n, p, alpha, Rninv);
	complex_comp_sig2(n, p, q, yR, yI, X, beta, theta, Rninv, &sig2);
	LL_old = compute_logL(p, n, alpha, sig2, 0);
	while(conv==0 && iter<max_iter){
		++iter;
		complex_res_beta_theta(n, q, p, m, C, X, yR, yI, beta, 
			&theta, 0, Rninv);
		complex_comp_alpha(n, p, D, alpha);
		comp_Rninv(n, p, alpha, Rninv);
		complex_comp_sig2(n, p, q, yR, yI, X, beta, theta, Rninv, &sig2);
		LL_new = compute_logL(p, n, alpha, sig2, 0);
		if(LL_new - LL_old < LL_eps) conv =1;
		else LL_old = LL_new;
	}
	store_est(q, p, beta, theta, sig2, alpha, par);
	FREE_MATRIX(D);
	FREE_MATRIX(Rninv);
	FREE_VECTOR(beta);
	FREE_VECTOR(alpha);	
}

void complex_res_beta_theta(int n, int q, int p, int m, double **C, 
	double **X, double *yR, double *yI, double *beta_til, 
	double *theta_til, int Rn_is_In, double **Rninv)
{
	int i;
	double num, den;
	double *bR;		MAKE_VECTOR(bR, q);
	double *bI;		MAKE_VECTOR(bI, q);
	double **Psi;		MAKE_MATRIX(Psi, q, q);
	double *Psi_bR;	MAKE_VECTOR(Psi_bR, q);
	double *Psi_bI;	MAKE_VECTOR(Psi_bI, q);
	double **XRX;		MAKE_MATRIX(XRX, q, q);
	GLS(n, q, p, X, Rninv, yR, bR);
	GLS(n, q, p, X, Rninv, yI, bI);
	comp_Psi(n, q, p, m, X, Rninv, Rn_is_In, C, Psi);
	matxvec(Psi, q, q, bR, q, Psi_bR);
	matxvec(Psi, q, q, bI, q, Psi_bI);
	comp_XRX(n, q, p, X, Rninv, XRX);
	num = 2.0 * bilinear(Psi_bR, XRX, bI, q);
	den = bilinear(Psi_bR, XRX, bR, q) - bilinear(Psi_bI, XRX, bI, q);
	*theta_til = 0.5 * atan2(num, den);
	for(i=0; i<q; ++i)	
		beta_til[i] = Psi_bR[i]*cos(*theta_til) + Psi_bI[i]*sin(*theta_til);
	if(beta_til[0]<0){
		for(i=0; i<q; ++i)	beta_til[i] = -1.*beta_til[i];
		if(*theta_til>0)		
			*theta_til = *theta_til - M_PI;
		else
			*theta_til = *theta_til +M_PI;
	}
	FREE_VECTOR(bR);
	FREE_VECTOR(bI);
	FREE_MATRIX(Psi);
	FREE_VECTOR(Psi_bR);
	FREE_VECTOR(Psi_bI);
	FREE_MATRIX(XRX);
}

void compute_lrt_stat(int n, int p, int q, int complex_mag, double *par_unres, 
	double *par_res, double *lrt_stat)
{
	int i;
	double *alpha;		MAKE_VECTOR(alpha, p);
	double sig2, logL_unres, logL_res;
	//logL_unres
	sig2 = par_unres[q+1];
	for(i=0; i<p; ++i)	alpha[i] = par_unres[q+2+i];
	logL_unres = compute_logL(p, n, alpha, sig2, complex_mag);
	//logL_res
	sig2 = par_res[q+1];
	for(i=0; i<p; ++i)	alpha[i] = par_res[q+2+i];
	logL_res = compute_logL(p, n, alpha, sig2, complex_mag);
	*lrt_stat = 2.0 * (logL_unres - logL_res);
	FREE_VECTOR(alpha);
}

void store_est(int q, int p, double *beta, double theta, double sig2, 
	double *alpha, double *par)
{
	int i;
	for(i=0; i<q; ++i)	par[i] = beta[i];
	par[q] = theta;
	par[q+1] = sig2;
	for(i=0; i<p; ++i)	par[q+2+i] = alpha[i];
}

// For order detection

void Rwrapper_complex_unres_only(int *n, int *q, int *p, double *X_vec, 
	double *yR, double *yI, int *max_iter, double *LL_eps, 
	double *par, double *LL)
{
	double **X;			MAKE_MATRIX(X, *n, *q);
	copy_1d_to_2d(X_vec, X, *n, *q);
	complex_unres_only(*n, *q, *p, X, yR, yI, *max_iter, *LL_eps, 
		par, LL);
	FREE_MATRIX(X);
}

void complex_unres_only(int n, int q, int p, double **X, 
	double *yR, double *yI, int max_iter, double LL_eps, 
	double *par, double *LL)
{
	int i;
	double *beta;		MAKE_VECTOR(beta, q);
	double *alpha;	MAKE_VECTOR(alpha, p);
	double theta, sig2;
	if(p==0){
		complex_indep_unres(X, yR, yI, n, q, beta, &sig2, &theta);
		for(i=0; i<q; ++i) par[i] = beta[i];
		par[q] = theta;
		par[q+1] = sig2;
		*LL = -1.0 * n * log(sig2) - n;
	}
	else{
		complex_unres_est(n, q, p, X, yR, yI, max_iter, par, 
			LL_eps);
		sig2 = par[q+1];
		for(i=0; i<p; ++i)	alpha[i] = par[q+2+i];
		*LL = compute_logL(p, n, alpha, sig2, 0);
	}
	FREE_VECTOR(beta);
	FREE_VECTOR(alpha);
}