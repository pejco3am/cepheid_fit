#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cmpfit-1.2/mpfit.h"
#include "constants.h"
#include "data_structure.h"
#include "data_managment.h"


void initialize_data_structure(data_structure *data) {
	data->id = new int[NDATA_MAX];
	data->iflt = new int[NDATA_MAX];
	data->iref = new int[NDATA_MAX];
	data->jd = new double[NDATA_MAX];
	data->mag = new double[NDATA_MAX];
	data->err = new double[NDATA_MAX];
	
	data->nperflt = new int[NFLT_MAX];
	for (int i=0;i<NFLT_MAX;i++) data->nperflt[i] = 0;
}



void initialize_fitting_coefs(cepheid *cep, double *a) {
	// the order of coefficients is period pars + amplitude + mbar + beta + fourier
	for (int i=0;i<N_PER_PARS_MAX;i++) a[i] = cep->e.per_pars[i];
	a[index_amp()] = cep->p.a;
	a[index_rho0()] = cep->p.rho0;
	for (int i=0;i<NFLT_MAX;i++) {
		a[index_mbar() + i] = cep->p.mbar[i];
		a[index_beta() + i] = cep->p.beta[i];
	}
	for (int i=0;i<(4*NF_MAX);i++) {
		a[index_fc() + i] = cep->p.fc[i];
	}
}


void initialize_cep_pars(cep_pars *cep) {
	cep->nf = NF_MAX;
	cep->fc = new double[4*cep->nf];
	for (int i=0;i<(4*cep->nf);i++) cep->fc[i] = 0.0;
	cep->fc[0] = 1.0;	// normalization of the template
	cep->fc[NF_MAX] = 0.0;
	cep->a = 0.24;
	cep->rho0 = 0.1;
	cep->mbar = new double[NFLT_MAX];
	cep->beta = new double[NFLT_MAX];

	for (int i=0;i<NFLT_MAX;i++) {
		cep->mbar[i] = 0.0;
		cep->beta[i] = 1.0;
	}


}


void initialize_ephemeris(ephemeris *eph) {
		for (int i=0;i<NCEP_MAX;i++) {
			eph[i].mode = -10;		// negative value indicates info not available
			eph[i].per_pars = new double[N_PER_PARS_MAX];
			eph[i].per_pars[0] = -10.00;	// negative period indicates info not available
			for (int j=1;j<N_PER_PARS_MAX;j++) eph[i].per_pars[j] = 0.0;	// zero all other parameters
	}
}

void initialize_fitting_pars_single(mp_par *pars, const int ncoef) {
	// first set all pars to be fixed
	for (int i=0;i<ncoef;i++) pars[i].fixed = 1;
}

void vary_phase(mp_par *pars) {
	pars[1].fixed = 0;
}

void vary_amp(mp_par *pars) {
	pars[index_amp()].fixed = 0;
}

void vary_rho0(mp_par *pars, cepheid *cep) {
	if (cep->d.nperflt[0] > 0) {
		pars[index_rho0()].fixed = 0;
	}
}

void vary_meanRV(mp_par *pars, cepheid *cep) {
	if (cep->d.nperflt[0] > 0) {
		pars[index_mbar()+0].fixed = 0;
	}
}

void vary_mbar(mp_par *pars, cepheid *cep) {
	// switch on mbar and beta only if there are actually data to constrain
	for (int i=0;i<NFLT_MAX;i++) if (cep->d.nperflt[i] > 0) {
			pars[index_mbar() + i].fixed = 0;
	}
}

void vary_beta_specific(mp_par *pars, const int iflt) {
	pars[index_beta()+iflt].fixed = 0;
}

void vary_fc(mp_par *pars, cepheid *cep) {
	// now vary all coefficients of the Fourier expansion up to order NF, but don't change the first one, which sets the overal normalization
	for (int i=0;i<cep->p.nf;i++) {
		pars[index_fc() + i].fixed = 0;
		pars[index_fc() + NF_MAX + i].fixed = 0;
		pars[index_fc() + 2*NF_MAX + i].fixed = 0;
		pars[index_fc() + 3*NF_MAX + i].fixed = 0;	
	}
	pars[index_fc()].fixed = 1;		// fix the zeroth order coeffs for radius
	pars[index_fc()+NF_MAX].fixed = 1;

}
