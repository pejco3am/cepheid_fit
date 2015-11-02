void initialize_data_structure(data_structure *data);
int max_fitted_pars(data_structure *data);
void initialize_ephemeris(ephemeris *eph);
void initialize_cep_pars(cep_pars *cep);
void initialize_fitting_coefs(cepheid *cep, double *a);
void initialize_fitting_pars_single(mp_par *pars, const int ncoef);
void vary_phase(mp_par *pars);
void vary_amp(mp_par *pars);
void vary_rho0(mp_par *pars, cepheid *cep);
void vary_meanRV(mp_par *pars, cepheid *cep);
void vary_mbar(mp_par *par, cepheid *cep);
void vary_beta_specific(mp_par *par, const int iflt);
void vary_fc(mp_par *pars, cepheid *cep);
