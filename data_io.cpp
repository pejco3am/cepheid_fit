#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "data_structure.h"

void read_data(data_structure *data, char const  *filname) {
	int i = 0;
	FILE *dov = fopen(filname, "r");
	char dummy[1000];
	while (fgets(dummy, 1000, dov)) {
		sscanf(dummy, "%i %i %lf  %lf %lf %i", &data->iflt[i], &data->iref[i], &data->jd[i], &data->mag[i], &data->err[i], &data->id[i]);

		// put in some reasonable uncertainties, if missing
		if (data->iflt[i] == 0)  data->err[i] = fmax(data->err[i], 1.0);	// 1km/s
		if (data->iflt[i] > 0) data->err[i] = fmax(data->err[i], 0.02);		// 0.02 mag
		data->nperflt[data->iflt[i]]++;
		i++;
	}
	fclose(dov);
	data->ndata = i;
}

void read_periods(ephemeris *eph, char const *filname) {
	int id, mode;
	double a[10];
	FILE *dov = fopen(filname, "r");
	char dummy[1000], dummy2[100];
	while (fgets(dummy, 1000, dov)) {
		// hard-code 10 period pars per cepheid
		sscanf(dummy, "%i %s %i  %le %le %le %le %le %le %le %le %le %le", &id, &dummy2, &mode, &a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7], &a[8], &a[9]);
		eph[id].mode = mode;
		strcpy(eph[id].name, dummy2);
		// copy cepheid period pars
		for (int i=0;i<10;i++) eph[id].per_pars[i] = a[i];
	}
	fclose(dov);

}

void write_ephemeris_line(ephemeris *eph, FILE *ven, const int id) {
	fprintf(ven, "%i %s %i  %.7f", id, eph->name, eph->mode, eph->per_pars[0]);
	for (int i=1;i<N_PER_PARS_MAX;i++) fprintf(ven, " %.3le", eph->per_pars[i]);
	fprintf(ven, "\n");
}


void write_periods(ephemeris *eph, char const *filname, const int id) {
	// only modify one record
	
	// open a temporary file 
	FILE *temp = fopen("gafsjglfk.dat", "w");
	FILE *dov = fopen(filname, "r");

	char dummy[10000];
	int cep_id, written = 0;
	while (fgets(dummy, 10000, dov)) {
		sscanf(dummy, "%i", &cep_id);
		if (cep_id != id) {
			fputs(dummy, temp);
		} else {	
			// now write the results in the file
			write_ephemeris_line(eph, temp, id);
			written = 1;
		}
	}
	if (written == 0) write_ephemeris_line(eph, temp, id);
	fclose(dov);
	fclose(temp);
	remove(filname);
	rename("gafsjglfk.dat", filname);
}


void read_beta(double *beta, char const *filname) {
	FILE *dov = fopen(filname, "r");
	int i =1;
	char dummy[1000];

	while (fgets(dummy, 1000, dov)) {
			// hard-code 10 period pars per cepheid
		sscanf(dummy, "%lf", &beta[i]);
		i++;
	}
}

void read_fc(double *fc, char const *filname, const double period) {
	double lper = log10(period/10.0);


	FILE *dov = fopen(filname, "r");
	char dummy[10000], *dummy2;
	int NF, NP, offset, lper_index, dummy_int;
	fgets(dummy, 1000, dov);
	sscanf(dummy, "%i %i", &NF, &NP);
	double *lperbins = new double[NP];	

	// now read in the period bins
	fgets(dummy, 1000, dov);
	dummy2 = dummy;
	for (int i=0;i<NP;i++) {
		sscanf(dummy2, " %lf%n", &lperbins[i], &offset);
		if (lper > lperbins[i]) lper_index = i;
		dummy2 += offset;
	}
	// setup interpolation coefficients;
	double w, ct, st;
	if (lper < lperbins[0]) {
		lper_index = 0;
		w = 0;
	} else if (lper > lperbins[NP-1]) {
		lper_index = NP-2;
		w = 1;
	} else {
		w = (lper-lperbins[lper_index])/(lperbins[lper_index+1]-lperbins[lper_index]);
	}

	printf("%f %f %i\n", lper, w, lper_index);

	// read in temperature coefs
	for (int iline=0;iline<NF;iline++) {
		fgets(dummy, 10000, dov);
		dummy2 = dummy;
		sscanf(dummy2, " %i%n", &dummy_int, &offset);
		dummy2 += offset;

		fc[2*NF_MAX+iline] = 0;
		fc[3*NF_MAX+iline] = 0;

		for (int i=0;i<NP;i++) {
			sscanf(dummy2, " %le %le%n", &ct, &st, &offset);
			if (i == lper_index) {
				fc[2*NF_MAX+iline] += (1-w)*ct;
				fc[3*NF_MAX+iline] += (1-w)*st;
			}
			if (i == (lper_index+1)) {
				fc[2*NF_MAX+iline] += w*ct;
				fc[3*NF_MAX+iline] += w*st;
			}
			dummy2 += offset;
		}
	}

	// read in radius coefs
	for (int iline=0;iline<NF;iline++) {
		fgets(dummy, 10000, dov);
		dummy2 = dummy;
		sscanf(dummy2, " %i%n", &dummy_int, &offset);
		dummy2 += offset;

		fc[iline] = 0;
		fc[NF_MAX+iline] = 0;
		for (int i=0;i<NP;i++) {
			sscanf(dummy2, " %le %le%n", &ct, &st, &offset);
			if (i == lper_index) {
				fc[iline] += (1-w)*ct;
				fc[NF_MAX+iline] += (1-w)*st;
			}
			if (i == (lper_index+1)) {
				fc[iline] += w*ct;
				fc[NF_MAX+iline] += w*st;
			}
			dummy2 += offset;
		}
	}


	delete[] lperbins;
}

void write_cep_pars_line(cep_pars *cep, FILE *ven, const int id) {
	fprintf(ven, "%i %.3e %.3e %i %i", id, cep->a, cep->rho0, NFLT_MAX, NF_MAX);
	for (int i=0;i<NFLT_MAX;i++) fprintf(ven, " %.3e", cep->mbar[i]);
	for (int i=0;i<NFLT_MAX;i++) fprintf(ven, " %.3e", cep->beta[i]);
	for (int i=0;i<(4*NF_MAX);i++) fprintf(ven, " %.3e", cep->fc[i]);
	fprintf(ven, "\n");
}


void write_cep_pars(cep_pars *cep, char const *filname, const int id) {
	// only modify one record
	
	// open a temporary file 
	FILE *temp = fopen("gafsjglfk.dat", "w");
	FILE *dov = fopen(filname, "r");

	char dummy[10000];
	int cep_id, written = 0;
	while (fgets(dummy, 10000, dov)) {
		sscanf(dummy, "%i", &cep_id);
		if (cep_id != id) {
			fputs(dummy, temp);
		} else {	
			// now write the results in the file
			write_cep_pars_line(cep, temp, id);
			written = 1;
		}
	}
	if (written == 0) write_cep_pars_line(cep, temp, id);
	fclose(dov);
	fclose(temp);
	remove(filname);
	rename("gafsjglfk.dat", filname);
}

int read_cep_pars(cep_pars *cep, char const *filname, const int id) {
	FILE *dov = fopen(filname, "r");
	char dummy[10000];
	char *dummy2;
	int cep_id, offset, loc_NFLT_MAX, loc_NF_MAX;
	while (fgets(dummy, 10000, dov)) {
		dummy2 = dummy;
		sscanf(dummy2, "%i%n", &cep_id, &offset);
		if (cep_id != id) {continue;} else {
			dummy2 += offset;
	
			sscanf(dummy2, " %le%n", &cep->a, &offset);
			dummy2 += offset;

			sscanf(dummy2, " %le%n", &cep->rho0, &offset);
			dummy2 += offset;

			sscanf(dummy2, " %i %i%n", &loc_NFLT_MAX, &loc_NF_MAX, &offset);
			dummy2 += offset;


			if (loc_NFLT_MAX > NFLT_MAX) {
				printf("increase NFLT_MAX in constants.h\n");
				exit(0);
			}
			if (loc_NF_MAX > NF_MAX) {
				printf("increase NF_MAX in constants.h\n");
				exit(0);
			}

			for (int i=0;i<loc_NFLT_MAX;i++) {
				sscanf(dummy2, " %le%n", &cep->mbar[i], &offset);
				dummy2 += offset;
			}
			
			for (int i=0;i<loc_NFLT_MAX;i++) {
				sscanf(dummy2, " %le%n", &cep->beta[i], &offset);
				dummy2 += offset;
			}


			for (int i=0;i<(4*loc_NF_MAX);i++) {
				sscanf(dummy2, " %le%n", &cep->fc[i], &offset);
				dummy2 += offset;
			}

			return 1;
		}
	}
	return 0;
}


