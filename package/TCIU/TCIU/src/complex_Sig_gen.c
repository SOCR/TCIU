#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "array.h"
#include "mat_vec.h"
#include "RC_interface.h"
#include "GLS.h"
#include "complex_Sig_gen.h"

void Rwrapper_est_complex_par_ri_temp_dep(int *n, int *q, int *p, double *yr, 
	double *yi, double *X_vec, double *beta, double *theta, double *sr2, 
	double *si2, double *rho, double *alpha, double *tol, int *max_iter, 
	double *LL_opt)
{
	
	double **X;			MAKE_MATRIX(X, *n, *q);
	copy_1d_to_2d(X_vec, X, *n, *q);
	est_complex_par_ri_temp_dep(*n, *q, *p, yr, yi, X, beta, theta, 
		sr2, si2, rho, alpha, *tol, *max_iter, LL_opt);
	FREE_MATRIX(X);
}

void complex_Sig_gen_p0(int n, int q, double *yr, double *yi, double **X, 
	double *beta, double *theta, double *sr2, double *si2, double *rho, 
	double *LL_opt)
{
	compute_beta_theta_indep(n, q, yr, yi, X, beta, theta);
	comp_Sigma_indep(n, q, yr, yi, X, beta, *theta, sr2, si2, rho);
	
}

void est_complex_par_ri_temp_dep(int n, int q, int p, double *yr, 
	double *yi, double **X, double *beta, double *theta, double *sr2, 
	double *si2, double *rho, double *alpha, double tol, int max_iter, 
	double *LL_opt)
{
	int conv=0, iter=0;
	double LL_old, LL_new;
	compute_beta_theta_indep(n, q, yr, yi, X, beta, theta);
	comp_Sigma_indep(n, q, yr, yi, X, beta, *theta, sr2, si2, rho);
	complex_comp_alpha_ri_dep(n, p, q, X, beta, *theta, yr, yi, 
		*sr2, *si2, *rho, alpha);
	compute_LL_ri_time_dep(n, p, *sr2, *si2, *rho, alpha, &LL_new);
	while(conv == 0 && iter < max_iter){
		LL_old = LL_new;
		++iter;
		update_beta_theta_Sigma(n, q, p, yr, yi, X, beta, theta, sr2, 
			si2, rho, alpha);
		complex_comp_alpha_ri_dep(n, p, q, X, beta, *theta, yr, yi, 
			*sr2, *si2, *rho, alpha);
		compute_LL_ri_time_dep(n, p, *sr2, *si2, *rho, alpha, 
			&LL_new);
		if(LL_new - LL_old < tol) conv = 1;
	}
	*LL_opt = LL_new;
}

void compute_LL_ri_time_dep(int n, int p, double sr2, 
	double si2, double rho, double *alpha, double *LL)
{
	double **Rpinv;		MAKE_MATRIX(Rpinv, p, p);
	comp_Rpinv(p, alpha, Rpinv); //GLS.c
	double det = fabs(my_det(p, Rpinv));
	*LL = -1.0* n / 2 * log(sr2 * si2 * (1 - rho*rho)) + log(det) - n;
	FREE_MATRIX(Rpinv);
}

void compute_beta_theta_indep(int n, int q, double *yr, double *yi, double **X, 
	double *beta, double *theta)
{
	int i;
	double num, den;
	double *bR;		MAKE_VECTOR(bR, q);
	double *bI;		MAKE_VECTOR(bI, q);
	double **XpX;		MAKE_MATRIX(XpX, q, q); 
	OLS(X, yr, n, q, bR);
	OLS(X, yi, n, q, bI);
	X_transpose_X(X, n, q, XpX);
	num = 2.*bilinear(bR, XpX, bI, q);
	den = quadratic(XpX, bR, q) - quadratic(XpX, bI, q);
	*theta = 0.5*atan2(num, den);
	for(i=0; i<q; ++i)	beta[i] = bR[i]*cos(*theta)+bI[i]*sin(*theta);
	FREE_VECTOR(bR);
	FREE_VECTOR(bI);
	FREE_MATRIX(XpX);
}

void comp_Sigma_indep(int n, int q, double *yr, double *yi, double **X, 
	double *beta, double theta, double *sr2, double *si2, double *rho)
{
	int t;
	double temp=0.0;
	double *mu;		MAKE_VECTOR(mu, n);
	matxvec(X, n, q, beta, q, mu);
	for(t=0; t<n; ++t) 	temp += pow(yr[t] - mu[t] * cos(theta), 2);
	*sr2 = temp / n; temp=0.0;
	for(t=0; t<n; ++t) 	temp += pow(yi[t] - mu[t] * sin(theta), 2);
	*si2 = temp / n; temp=0.0;
	for(t=0; t<n; ++t)
		temp += (yr[t] - mu[t] * cos(theta)) * (yi[t] - mu[t] * sin(theta));
	*rho = temp / (n * sqrt(*sr2 * *si2));
	FREE_VECTOR(mu);
}



void complex_comp_alpha_ri_dep(int n, int p, int q, double **X, 
	double *beta, double theta, double *yr, double *yi, double sr2, 
	double si2, double rho, double *alpha)
{
	int i, j;
	double **D;			MAKE_MATRIX(D, p+1, p+1);
	double **alp_mat;	MAKE_MATRIX(alp_mat, p, p);
	double *d;			MAKE_VECTOR(d, p);
	complex_compD_ri_dep(n, p, q, X, beta, theta, yr, yi, sr2, si2, 
		rho, D);
	//solve for alpha
	for(i=1; i<=p; ++i)	d[i-1] = D[i][0];
	for(i=1; i<=p; ++i)
		for(j=1; j<=p; ++j)		
			alp_mat[i-1][j-1] = D[i][j] + j/n * D[0][abs(j-i)];
	my_inv(p, alp_mat);
	matxvec(alp_mat, p, p, d, p, alpha);
	FREE_MATRIX(D);
	FREE_MATRIX(alp_mat);
	FREE_VECTOR(d);
}

void complex_compD_ri_dep(int n, int p, int q, double **X, double *beta, 
	double theta, double *yr, double *yi, double sr2, double si2, 
	double rho, double **D)
{
	int t, i, k;
	double drr = 0.0, dii=0.0, dri=0.0, dir = 0.0, 
		cosine=cos(theta), sine=sin(theta);
	double *mu;		MAKE_VECTOR(mu, n);
	double *er;		MAKE_VECTOR(er, n);
	double *ei;		MAKE_VECTOR(ei, n);
	matxvec(X, n, q, beta, q, mu);
	for(t=0; t<n; ++t){
		er[t] = yr[t] - mu[t] * cosine;
		ei[t] = yi[t] - mu[t] * sine;
	}
	for(i=0; i<=p; ++i){
		for(k=0; k<=p-i; ++k){
			for(t=i; t<n-i-k; ++t){
				drr += er[t] * er[t+k];
				dii += ei[t] * ei[t+k];
				dri += er[t] * ei[t+k];
				dir += ei[t] * er[t+k];
			}
			D[i][i+k] = D[i+k][i] = drr / sr2 + dii / si2 - 
									rho / sqrt(sr2*si2) * (dri + dir);
			drr = dii = dri = dir = 0.0;
		}
	}
	FREE_VECTOR(mu);
	FREE_VECTOR(er);
	FREE_VECTOR(ei);

	
}

void update_beta_theta_Sigma(int n, int q, int p, double *yr, 
	double *yi, double **X, double *beta, double *theta, double *sr2, 
	double *si2, double *rho, double *alpha)
{
	double *br;		MAKE_VECTOR(br, q);
	double *bi;		MAKE_VECTOR(bi, q);
	double **Rninv;	MAKE_MATRIX(Rninv, n, n);
	comp_Rninv(n, p, alpha, Rninv); //GLS.c
	GLS(n, q, p, X, Rninv, yr, br); //GLS.c
	GLS(n, q, p, X, Rninv, yi, bi);
	update_theta(n, q, p, yr, yi, X, br, bi, theta, *sr2, *si2, 
		*rho, Rninv);
	update_beta(n, q, p, yr, yi, X, br, bi, *theta, *sr2, *si2, 
		*rho, beta);
	update_Sigma(n, q, p, yr, yi, X, beta, *theta, sr2, si2, rho, 
		Rninv);
	FREE_VECTOR(br);
	FREE_VECTOR(bi);
	FREE_MATRIX(Rninv);
}

void update_theta(int n, int q, int p, double *yr, 
	double *yi, double **X, double *br, double *bi, double *theta, double sr2, 
	double si2, double rho, double **Rninv)
{
	double **XRX;	MAKE_MATRIX(XRX, q, q);
	comp_XRX(n, q, p, X, Rninv, XRX); //GLS.c
	double Brr = quadratic(XRX, br, q);
	double Bii = quadratic(XRX, bi, q);
	double Bri = bilinear(br, XRX, bi, q);
	double D = Brr/(sr2*sr2) + rho*rho/(sr2*si2)*Bii - 
		2*rho/(pow(sr2,1.5) * sqrt(si2)) * Bri;
    double E = Bii/(si2*si2) + rho*rho/(sr2*si2)*Brr - 
		2*rho/(sqrt(sr2) * pow(si2,1.5)) * Bri;
    double F = Bri*(1+rho*rho)/(sr2*si2) - rho/sqrt(sr2*si2) * (Brr/sr2 + Bii/si2);
    double A = D/si2 - E/sr2;
    double B = -1.0*F * (1/sr2 + 1/si2) - rho/sqrt(sr2*si2) * (D+E);
    double C = rho/sqrt(sr2*si2) * (D-E) + F*(1/sr2 - 1/si2);
	double psi = atan2(B, A);
    *theta = 0.5 * (asin(C/sqrt(A*A+B*B)) - psi);
	FREE_MATRIX(XRX);
}

void update_beta(int n, int q, int p, double *yr, 
	double *yi, double **X, double *br, double *bi, double theta, 
	double sr2, double si2, double rho, double *beta)
{
	int i;
	double numr = cos(theta) / sr2 - rho * sin(theta) / sqrt(sr2*si2);
	double numi = sin(theta) / si2 - rho * cos(theta) / sqrt(sr2*si2);
	double den = cos(theta)*cos(theta)/sr2 + sin(theta)*sin(theta)/si2 - 
				2.0*rho/sqrt(sr2*si2) * sin(theta) * cos(theta);
	for(i=0; i<q; ++i) 
		beta[i] = (br[i] * numr + bi[i] * numi) / den;
}

void update_Sigma(int n, int q, int p, double *yr, double *yi, 
	double **X, double *beta, double theta, double *sr2, 
	double *si2, double *rho, double **Rninv)
{
	int t, tp;
	double *mu;		MAKE_VECTOR(mu, n);
	double *er;		MAKE_VECTOR(er, n);
	double *ei;		MAKE_VECTOR(ei, n);
	double cosine=cos(theta), sine=sin(theta);
	double sumrr=0.0, sumii=0.0, sumri=0.0;
	matxvec(X, n, q, beta, q, mu);
	for(t=0; t<n; ++t){
		er[t] = yr[t] - mu[t] * cosine;
		ei[t] = yi[t] - mu[t] * sine;
	}
	for(t=0; t<n; ++t){
		for(tp=max(t-p,0); tp<=min(t+p,n-1); ++tp){
			sumrr += er[t]*Rninv[t][tp]*er[tp];
			sumii += ei[t]*Rninv[t][tp]*ei[tp];
			sumri += er[t]*Rninv[t][tp]*ei[tp];
		}
	}
	*sr2 = sumrr / n;
	*si2 = sumii / n;
	*rho = sumri / (n * sqrt(*sr2 * *si2));
	FREE_VECTOR(mu);
	FREE_VECTOR(er);
	FREE_VECTOR(ei);
}

