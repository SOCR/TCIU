void Rwrapper_complex_activation(int *n, int *q, int *p, int *m, 
	double *C_vec, double *X_vec, double *yR, double *yI, 
	int *max_iter, double *LL_eps, double *par_vec, 
	double *lrt_stat);

void complex_activation(int n, int q, int p, int m, double **C, 
	double **X, double *yR, double *yI, int max_iter, 
	double LL_eps, double *par_unres, double *par_res, double *lrt_stat);

//Independent functions

void complex_indep(int n, int q, int m, double **C, double **X, 
	double *yR, double *yI, double *par_unres, double *par_res,
	double *lrt_stat);

void complex_indep_unres(double **X, double *yR, double *yI, int n, int q, 
	double *beta, double *sig2, double *theta);

void complex_res_indep(int n, int q, int m, double **C, double **X, 
	double *yR, double *yI, double *beta, double *sig2, double *theta);

void comp_Psi(int n, int q, int p, int m, double **X, double **Rninv, 
	int Rn_is_In, double **C, double **Psi);

// AR(p) functions

void complex_arp(int n, int q, int p, int m, double **X, double **C, 
	double *yR, double *yI, int max_iter, double LL_eps, 
	double *par_unres, double *par_res, double *lrt_stat);

void complex_unres_est(int n, int q, int p, double **X, 
	double *yR, double *yI, int max_iter,
	double *par, double LL_eps);

void comp_theta_beta_arp(int n, int q, int p, double **X, 
	double *yr, double *yi, double **Rninv,
	double *beta, double *theta);

void complex_compD(int n, int p, int q, double **X, double *beta, double theta, double *yr,
	double *yi, double **D);

void complex_comp_alpha(int n, int p, double **D, double *alpha);

void complex_comp_sig2(int n, int p, int q, double *yR, double *yI, double **X, double *beta, 
	double theta, double **Rninv, double *sig2);

double compute_logL(int p, int n, double *alpha, double sig2, 
	int complex_mag);

void complex_res_est(int n, int q, int p, int m, double **X, 
	double **C, double *yR, double *yI, int max_iter,
	double *par, double LL_eps);

void complex_res_beta_theta(int n, int q, int p, int m, double **C, 
	double **X, double *yR, double *yI, double *beta_til, 
	double *theta_til, int Rn_is_In, double **Rninv);

void compute_lrt_stat(int n, int p, int q, int complex_mag, double *par_unres, 
	double *par_res, double *lrt_stat);

void store_est(int q, int p, double *beta, double theta, double sig2, 
	double *alpha, double *par);

// For order detection

void Rwrapper_complex_unres_only(int *n, int *q, int *p, double *X_vec, 
	double *yR, double *yI, int *max_iter, double *LL_eps, 
	double *par, double *LL);

void complex_unres_only(int n, int q, int p, double **X, 
	double *yR, double *yI, int max_iter, double LL_eps, 
	double *par, double *LL);
	