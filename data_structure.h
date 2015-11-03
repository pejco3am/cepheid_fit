struct ephemeris {
	// information about period and phase
	int mode;	// to allow for different functional prescription of P(t)
	int n_per_pars;
	double *per_pars;
	char name[100];
};



struct cep_pars {
	// information about Fourier series
	int nf;	// order of Fourier series
	double *fc;

	double a;	// amplitude
	double rho0;	// mean radius
	double *mbar;	// mean magnitudes of filters (0=RV, 1=U, etc)
	double *beta;	// betas for individual photometric bands (0 undefined, 1=U, etc)
};


struct data_structure {
	int ndata;
	// each of these arrays has ndata elements
	int *id;	// cepheid id
	int *iflt;	// filter of each observations (0 = RV, 1=U, 2=B, etc)
	int *iref;	// id of data reference (to allow for possible shifts etc in the future);
	double *jd;	// HJD
	double *mag;	// measurements (magnitude or RV)
	double *err;	// uncertainty


	int *nperflt;	// number of measurements in a given band
};



struct cepheid {
	ephemeris e;
	cep_pars p;
	data_structure d;
};


