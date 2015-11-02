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

	data_structure all_data;
	initialize_data_structure(&all_data);

	read_data(&all_data, "all.dat");
	ephemeris *eph = new ephemeris[NCEP_MAX];
	initialize_ephemeris(eph);
	read_periods(eph, "period.dat");

	int cep_fit = 577;

	// initialize data
	data_structure data;
	initialize_data_structure(&data);
	select_data(&all_data, &data, cep_fit);

	printf("%f\n", eph[cep_fit].per_pars[0]);

	// initialize cepheid parameters
	cep_pars my_cep_pars;
	initialize_cep_pars(&my_cep_pars);
	int in_catalog = read_cep_pars(&my_cep_pars, "cep_pars.out", cep_fit);
	if (in_catalog == 0) {
		read_beta(my_cep_pars.beta, "beta.dat");
		my_cep_pars.beta[28] = 4.23;
		read_fc(my_cep_pars.fc, "template.dat", eph[cep_fit].per_pars[0]);
	}


	// maximum fourier order is 5
	my_cep_pars.nf = 3;


	cepheid my_cepheid;
	my_cepheid.p = my_cep_pars;
	my_cepheid.d = data;
	my_cepheid.e = eph[cep_fit];

	//printf("%i %f\n", my_cepheid.e.mode, my_cepheid.e.per_pars[0]);


	int ncoef = N_PER_PARS_MAX + 2 + 2*NFLT_MAX + 4*NF_MAX;		// N_PER_PARS_MAX + amplitude + mean_radius +  mbar + beta + fourier
	//printf("%i\n", ncoef);
	double *a = new double[ncoef];
	initialize_fitting_coefs(&my_cepheid, a);


	mp_par pars[ncoef];
	memset(&pars[0], 0, sizeof(pars));

	initialize_fitting_pars_single(pars, ncoef);
	vary_phase(pars);
	vary_amp(pars);
	vary_rho0(pars, &my_cepheid);
	vary_meanRV(pars, &my_cepheid);
	vary_mbar(pars, &my_cepheid);
	vary_fc(pars, &my_cepheid);





	mp_result result;
	memset(&result, 0, sizeof(result));
	mp_config config;
	memset(&config, 0, sizeof(config));
	config.maxiter=200;
	result.xerror = new double[ncoef];
	//result.covar = covar;


	int status;
	status = mpfit(calculate_residuals, my_cepheid.d.ndata, ncoef, a, pars, &config, (void *) &my_cepheid, &result);


	// copy the results back to the original structure
	copy_fit_results_back(&my_cepheid, a);




	FILE *vel = fopen("velocity.dat", "w");
	FILE *pho = fopen("photometry.dat", "w");

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
	for (int i=0;i<N_PER_PARS_MAX;i++) {
		printf("per par %2i = %f +/- %f %i\n", i, a[i], result.xerror[i], pars[i].fixed);
	}
	printf("A = %f +/- %f  %i        rho0 = %f +/- %f  %i\n", a[index_amp()], result.xerror[index_amp()], pars[index_amp()].fixed, a[index_rho0()], result.xerror[index_rho0()], pars[index_rho0()].fixed);
	for (int i=0;i<NFLT_MAX;i++) {
		if (my_cepheid.d.nperflt[i] < 1) continue;
		printf("mbar %2i = %f +/- %f  %i     beta %2i = %f +/- %f  %i       %f  %f  %i\n", i, a[index_mbar()+i], result.xerror[index_mbar()+i], pars[index_mbar()+i].fixed, 
				i, a[index_beta()+i], result.xerror[index_beta()+i], pars[index_beta()+i].fixed,
				diffs[i], sqr_diffs[i], my_cepheid.d.nperflt[i]);
	}

	for (int i=0;i<NF_MAX;i++) {
		printf("%2i  %f +/- %f  %i  \t    %f +/- %f  %i      %f +/- %f  %i      %f +/- %f  %i \n", i, a[index_fc()+i], result.xerror[index_fc()+i], pars[index_fc()+i].fixed, 
					a[index_fc()+NF_MAX+i], result.xerror[index_fc()+NF_MAX+i], pars[index_fc()+NF_MAX+i].fixed, 
					a[index_fc()+2*NF_MAX+i], result.xerror[index_fc()+2*NF_MAX+i], pars[index_fc()+2*NF_MAX+i].fixed, 
					a[index_fc()+3*NF_MAX+i], result.xerror[index_fc()+3*NF_MAX+i], pars[index_fc()+3*NF_MAX+i].fixed);
	}


	write_cep_pars(&my_cepheid.p, "cep_pars.out", cep_fit);
	write_periods(&my_cepheid.e, "period.dat", cep_fit);


}
