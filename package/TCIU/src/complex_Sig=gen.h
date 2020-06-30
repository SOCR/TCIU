void Rwrapper_est_complex_par_ri_temp_dep(int *n, int *q, int *p, double *yr, 
	double *yi, double *X_vec, double *beta, double *theta, double *sr2, 
	double *si2, double *rho, double *alpha, double *tol, int *max_iter, 
	double *LL_opt);

void est_complex_par_ri_temp_dep(int n, int q, int p, double *yr, 
	double *yi, double **X, double *beta, double *theta, double *sr2, 
	double *si2, double *rho, double *alpha, double tol, int max_iter, 
	double *LL_opt);

void compute_LL_ri_time_dep(int n, int p, double sr2, 
	double si2, double rho, double *alpha, double *LL);

void compute_beta_theta_indep(int n, int q, double *yr, double *yi, double **X, 
	double *beta, double *theta);

void comp_Sigma_indep(int n, int q, double *yr, double *yi, double **X, 
	double *beta, double theta, double *sr2, double *si2, double *rho);

void complex_comp_alpha_ri_dep(int n, int p, int q, double **X, 
	double *beta, double theta, double *yr, double *yi, double sr2, 
	double si2, double rho, double *alpha);

void complex_compD_ri_dep(int n, int p, int q, double **X, double *beta, 
	double theta, double *yr, double *yi, double sr2, double si2, 
	double rho, double **D);

void update_beta_theta_Sigma(int n, int q, int p, double *yr, 
	double *yi, double **X, double *beta, double *theta, double *sr2, 
	double *si2, double *rho, double *alpha);

void update_theta(int n, int q, int p, double *yr, 
	double *yi, double **X, double *br, double *bi, double *theta, double sr2, 
	double si2, double rho, double **Rninv);

void update_beta(int n, int q, int p, double *yr, 
	double *yi, double **X, double *br, double *bi, double theta, 
	double sr2, double si2, double rho, double *beta);

void update_Sigma(int n, int q, int p, double *yr, double *yi, 
	double **X, double *beta, double theta, double *sr2, 
	double *si2, double *rho, double **Rninv);
