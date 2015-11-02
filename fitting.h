int calculate_residuals(int n_tot, int n_coef, double *a, double *dy, double **derivs, void *vars);
void get_template(double *a, const int nf, const double phase, double &del_rho, double &del_tau, double &del_drho);
double get_phase(const int mode, double *per_pars, const double jd);
double get_velocity(const double vbar, const double rho0, const double amp, const double del_rho, const double del_drho, const int mode, double *per_pars);
double get_magnitude(const double mbar, const double beta, const double amp, const double del_rho, const double del_tau);
