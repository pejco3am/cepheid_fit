#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "data_structure.h"


void select_data(data_structure *all_data, data_structure *data, const int cep_id) {
	int inew = 0;
	for (int i=0; i<all_data->ndata;i++) {
		if (all_data->id[i] == cep_id) {
			data->id[inew] = all_data->id[i];
			data->iflt[inew] = all_data->iflt[i];
			data->iref[inew] = all_data->iref[i];
			data->jd[inew] = all_data->jd[i];
			data->mag[inew] = all_data->mag[i];
			data->err[inew] = all_data->err[i];
			data->nperflt[data->iflt[inew]]++;
			inew++;
		}
	}
	data->ndata = inew;
}


int index_amp() {
	return N_PER_PARS_MAX;
}

int index_rho0() {
	return N_PER_PARS_MAX + 1;
}

int index_mbar() {
	return N_PER_PARS_MAX + 2;
}

int index_beta() {
	return N_PER_PARS_MAX + 2 + NFLT_MAX;
}

int index_fc() {
	return N_PER_PARS_MAX + 2 + 2*NFLT_MAX;
}


void copy_fit_results_back(cepheid *cep, double *a) {
	cep->p.a = a[index_amp()];
	cep->p.rho0 = a[index_rho0()];
	for (int i=0;i<NFLT_MAX;i++) {
		cep->p.mbar[i] = a[index_mbar()+i];
		cep->p.beta[i] = a[index_beta()+i];
	}
	for (int i=0;i<(4*NF_MAX);i++) {
		cep->p.fc[i] = a[index_fc()+i];
	}

	for (int i=0;i<N_PER_PARS_MAX;i++) cep->e.per_pars[i] = a[i];

}

