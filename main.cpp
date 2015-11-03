#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cmpfit-1.2/mpfit.h"
#include "constants.h"
#include "data_structure.h"
#include "initialize.h"
#include "data_io.h"
#include "fitting.h"
#include "data_managment.h"

inline double SQR(const double x) {
	return x*x;
}

int main () {

	// establish structure for data (see data_structure.h for details) and read all data
	data_structure all_data;
	initialize_data_structure(&all_data);
	read_data(&all_data, "all.dat");


	// choose which Cepheid to work on - the correspondence between the number and the name can be found in all.dat or catalog.dat
	int cep_fit = 21;


	// read ephemeris information for the desired cepheid from the database
	ephemeris eph;
	initialize_ephemeris(&eph);
	read_periods(&eph, "period.dat", cep_fit);

	// initialize data and select only data pertinent to the Cepheid of interest
	data_structure data;
	initialize_data_structure(&data);
	select_data(&all_data, &data, cep_fit);

	printf("Period: %f\n", eph.per_pars[0]);

	// initialize cepheid parameters
	cep_pars my_cep_pars;
	initialize_cep_pars(&my_cep_pars);
	// try to read in cepheid parameters from the database
	int in_catalog = read_cep_pars(&my_cep_pars, "cep_pars.out", cep_fit);
	// if Cepheid not in the database, read betas from PK12, fix the F620M beta to some appropriate value, and interpolate in the original PK12 templates to get an initial guess
	// sometimes the fitting gets lost on Fourier coefficients (crazy numbers and uncertainties for the temperature template), 
	// then force it to read PK12 templates by commenting the "if" part and do not vary Fourier coefficients (see below).
	// After A, rho0, etc is well established with fixed templates, you can vary the Fourier coefficients to improve the fit - this usually works
	if (in_catalog == 0) {
		read_beta(my_cep_pars.beta, "beta.dat");
		my_cep_pars.beta[28] = 4.23;
		read_fc(my_cep_pars.fc, "template.dat", eph.per_pars[0]);
	}


	// maximum fourier order
	my_cep_pars.nf = 10;

	cepheid my_cepheid;
	my_cepheid.p = my_cep_pars;
	my_cepheid.d = data;
	my_cepheid.e = eph;


	// Maximum number of coefficients to fit - see constants.h for definition

	int ncoef = N_PER_PARS_MAX + 2 + 2*NFLT_MAX + 4*NF_MAX;		// N_PER_PARS_MAX + amplitude + mean_radius +  mbar + beta + fourier
	double *a = new double[ncoef];
	initialize_fitting_coefs(&my_cepheid, a);

	mp_par pars[ncoef];
	memset(&pars[0], 0, sizeof(pars));

	// now pick which parts will we vary vary

	initialize_fitting_pars_single(pars, ncoef);		// everything fixed
	vary_phase(pars);					
	vary_amp(pars);						// sometimes useful to initially fix amplitude to some specific value
	vary_rho0(pars, &my_cepheid);
	vary_meanRV(pars, &my_cepheid);
	vary_mbar(pars, &my_cepheid);
	//vary_beta_specific(pars, 28);				// I found that varying beta for F620M doesn't the improve the fit dramatically
	vary_fc(pars, &my_cepheid);				// vary Fourier coefficients - useful to shutoff sometimes



	mp_result result;
	memset(&result, 0, sizeof(result));
	mp_config config;
	memset(&config, 0, sizeof(config));
	config.maxiter=200;
	result.xerror = new double[ncoef];

	int status;
	// actual fitting routine - see fitting.cpp for details
	status = mpfit(calculate_residuals, my_cepheid.d.ndata, ncoef, a, pars, &config, (void *) &my_cepheid, &result);


	// copy the results back to the original structure
	copy_fit_results_back(&my_cepheid, a);


	// now write the outputs
	FILE *vel = fopen("velocity.dat", "w");
	FILE *pho = fopen("photometry.dat", "w");


	// find standard deviations along the best-fit model
	double diffs[NFLT_MAX], sqr_diffs[NFLT_MAX];
	int nperflt[NFLT_MAX];
	for (int i=0;i<NFLT_MAX;i++) {
		diffs[i] = 0.0;
		sqr_diffs[i] = 0.0;
		nperflt[i] = 0;
	}


	for (int i=0;i<my_cepheid.d.ndata;i++) {
		double phase = get_phase(my_cepheid.e.mode, a, my_cepheid.d.jd[i]);
		double del_rho, del_tau, del_drho;
		get_template(a+index_fc(), my_cepheid.p.nf, phase, del_rho, del_tau, del_drho);
		double vmodel = get_velocity(a[index_mbar()], a[index_rho0()], a[index_amp()], del_rho, del_drho, my_cepheid.e.mode, a);
		int flt_id = my_cepheid.d.iflt[i];
		double mmodel = get_magnitude(a[index_mbar()+flt_id], a[index_beta()+flt_id], a[index_amp()], del_rho, del_tau);
		phase = phase - floor(phase);
		// print only data after JD 2 440 000
		if (my_cepheid.d.jd[i] < 40000) continue;
		if (my_cepheid.d.iflt[i]  == 0) {
			fprintf(vel, "%f %f %f %f %f   %.3e %.3e %.3e\n", my_cepheid.d.jd[i], phase, my_cepheid.d.mag[i], vmodel, my_cepheid.d.err[i], del_rho, del_drho, del_tau);
			diffs[0] += my_cepheid.d.mag[i]-vmodel;
			sqr_diffs[0] += SQR(my_cepheid.d.mag[i]-vmodel);
			nperflt[0]++;
		} else {
			fprintf(pho, "%i %f %f %f %f %f   %.3e %.3e\n", flt_id, my_cepheid.d.jd[i], phase, my_cepheid.d.mag[i], mmodel, my_cepheid.d.err[i], del_rho, del_tau);
			diffs[my_cepheid.d.iflt[i]] += my_cepheid.d.mag[i] - mmodel;
			sqr_diffs[my_cepheid.d.iflt[i]] += SQR(my_cepheid.d.mag[i] - mmodel);
			nperflt[my_cepheid.d.iflt[i]]++;
		}

	}
	fclose(vel);
	fclose(pho);

	// contains finely sampled model light curves
	FILE *mod = fopen("model.dat", "w");
	const int NPHASE = 200;
	for (int i=0;i<NPHASE;i++) {
		double phase = i/double(NPHASE);
		double del_rho, del_tau, del_drho;
		get_template(a+index_fc(), my_cepheid.p.nf, phase, del_rho, del_tau, del_drho);
		double vmodel = get_velocity(a[index_mbar()], a[index_rho0()], a[index_amp()], del_rho, del_drho, my_cepheid.e.mode, a);
		fprintf(mod, "%f  %f %f %f   %f", phase, del_rho, del_tau, del_drho,  vmodel);
		for (int i=1;i<NFLT_MAX;i++) {
			double mmodel = get_magnitude(a[index_mbar()+i], a[index_beta()+i], a[index_amp()], del_rho, del_tau);
			fprintf(mod, " %f", mmodel);
		}
		fprintf(mod, "\n");

	}
	fclose(mod);


	for (int i=0;i<NFLT_MAX;i++) {
		if (nperflt[i] > 1) {
			diffs[i] = diffs[i]/nperflt[i];
			sqr_diffs[i] = sqrt(  (sqr_diffs[i] - nperflt[i]*SQR(diffs[i]))/double(nperflt[i]-1.0));
		} else {
			diffs[i] = 0.0;
			sqr_diffs[i] = 0.0;
		}

	}




	printf("*** testlinfit status = %d, niter = %i\n", status, result.niter);
	printf("orig CHI2    = %f    (%d DOF)\n", result.orignorm, result.nfunc-result.nfree);
	printf("  CHI-SQUARE = %f    (%d DOF)\n", result.bestnorm, result.nfunc-result.nfree);
	printf("        NPAR = %d\n", result.npar);
	printf("       NFREE = %d\n", result.nfree);
	printf("     NPEGGED = %d\n", result.npegged);
	printf("     NITER = %d\n", result.niter);
	printf("      NFEV = %d\n", result.nfev);
	printf("       BIC = %f\n", result.bestnorm+ result.nfree*log(result.nfunc));
	printf("\n");
	printf("Period and phase parameters:\n");
	for (int i=0;i<my_cepheid.e.n_per_pars;i++) {
		printf("per par %2i = %f +/- %f %i\n", i, a[i], result.xerror[i], pars[i].fixed);
	}
	printf("Amplitude, rho0:\n");
	printf("A = %f +/- %f  %i        rho0 = %f +/- %f  %i\n", a[index_amp()], result.xerror[index_amp()], pars[index_amp()].fixed, a[index_rho0()], result.xerror[index_rho0()], pars[index_rho0()].fixed);
	printf("Parameters of individual bands (mbar, beta, and standard devitation along the best fit in each band:\n");
	for (int i=0;i<NFLT_MAX;i++) {
		if (my_cepheid.d.nperflt[i] < 1) continue;
		printf("mbar %2i = %f +/- %f  %i     beta %2i = %f +/- %f  %i     sigma  %i = %f  (N_obs=%i)\n", i, a[index_mbar()+i], result.xerror[index_mbar()+i], pars[index_mbar()+i].fixed, 
				i, a[index_beta()+i], result.xerror[index_beta()+i], pars[index_beta()+i].fixed,
				i, sqr_diffs[i], my_cepheid.d.nperflt[i]);
	}
	printf("Fourier parameters:\n");
	for (int i=0;i<my_cepheid.p.nf;i++) {
		printf("%2i  %f +/- %f  %i  \t    %f +/- %f  %i      %f +/- %f  %i      %f +/- %f  %i \n", i, a[index_fc()+i], result.xerror[index_fc()+i], pars[index_fc()+i].fixed, 
					a[index_fc()+NF_MAX+i], result.xerror[index_fc()+NF_MAX+i], pars[index_fc()+NF_MAX+i].fixed, 
					a[index_fc()+2*NF_MAX+i], result.xerror[index_fc()+2*NF_MAX+i], pars[index_fc()+2*NF_MAX+i].fixed, 
					a[index_fc()+3*NF_MAX+i], result.xerror[index_fc()+3*NF_MAX+i], pars[index_fc()+3*NF_MAX+i].fixed);
	}

	// write the fit results back to the database
	write_cep_pars(&my_cepheid.p, "cep_pars.out", cep_fit);
	write_periods(&my_cepheid.e, "period.dat", cep_fit);


}
