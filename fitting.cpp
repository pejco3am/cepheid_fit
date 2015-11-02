#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cmpfit-1.2/mpfit.h"
#include "constants.h"
#include "data_structure.h"
#include "data_managment.h"


double get_phase(const int mode, double *per_pars, const double jd) {
	double out_phase;
	switch (mode) {
		case 0:
			// simplest case of constant period
			out_phase = (jd-per_pars[1])/per_pars[0];
			break;
		case 1:
			// other possible period functional form
			break;

		default:
			out_phase = (jd-per_pars[1])/per_pars[0];
			break;
	}
	return out_phase;
}



void get_template(double *a, const int nf, const double phase, double &del_rho, double &del_tau, double &del_drho) {
	del_rho = 0;
	del_tau = 0;
	del_drho = 0;
	double ct = 1.0, st = 0.0, ctt, stt;
	double dct = cos(2*PI*phase), dst=sin(2*PI*phase);
	for (int i=0;i<nf;i++) {
		ctt = ct*dct - st*dst;
		stt = st*dct + ct*dst;
		//del_rho += a[i]*         cos(2*PI*(i+1)*phase) + 	a[NF_MAX+i]*  sin(2*PI*(i+1)*phase);
		//del_tau += a[2*NF_MAX+i]*cos(2*PI*(i+1)*phase) + 	a[3*NF_MAX+i]*sin(2*PI*(i+1)*phase);
		//del_drho+= 2*PI*(i+1)*(-a[i]*sin(2*PI*(i+1)*phase) + 	a[NF_MAX+i]*  cos(2*PI*(i+1)*phase));
		del_rho += a[i]*         ctt + 	a[NF_MAX+i]*  stt;
		del_tau += a[2*NF_MAX+i]*ctt + 	a[3*NF_MAX+i]*stt;
		del_drho+= 2*PI*(i+1)*(-a[i]*stt + 	a[NF_MAX+i]*ctt);
		ct = ctt;
		st = stt;
	}
}

double dphase_dt(const int mode, double *per_pars) {
	// should return derivative of phase with respect to time
	double out = 0.0;
	switch(mode) {
		case 0:
			// simplest case of constant period
			out = 1.0/per_pars[0];
			break;
		case 1:
			// some other period law
			break;
		default:
			out = 1.0/per_pars[0];
			break;
	}
	return out;
}


double get_velocity(const double vbar, const double rho0, const double amp, const double del_rho, const double del_drho, const int mode, double *per_pars) {
	// after period.f
//	double vcon = -R10day*ln10*amp/(projection_factor)*dphase_dt(mode, per_pars);
//	double vmod  = vcon*del_drho*exp(ln10*amp*del_rho);
//	return vbar + rho0*vmod;

	return vbar - rho0*R10day/projection_factor* exp(ln10*amp*amp*del_rho) * ln10*amp*amp*del_drho*dphase_dt(mode,per_pars);

}


double get_magnitude(const double mbar, const double beta, const double amp, const double del_rho, const double del_tau) {
	return mbar - 5.0*amp*amp*del_rho - 2.5*beta*amp*amp*del_tau;
}




// main fitting routine
int calculate_residuals(int n_tot, int n_coef, double *a, double *dy, double **derivs, void *vars) {
	struct cepheid *cep = (struct cepheid *) vars;

	for (int i=0;i<n_tot;i++) {
		dy[i] = 0.0;
		double phase = get_phase(cep->e.mode, a, cep->d.jd[i]);
		double del_rho, del_tau, del_drho;
		get_template(a+index_fc(), cep->p.nf, phase, del_rho, del_tau, del_drho);
		int flt_id = cep->d.iflt[i];
		//if (flt_id == 1) continue;
		if (cep->d.jd[i] < 40000) continue;
		if (flt_id == 0) {
			double vmodel = get_velocity(a[index_mbar()], a[index_rho0()], a[index_amp()], del_rho, del_drho, cep->e.mode, a);
			dy[i] = (cep->d.mag[i] - vmodel)/cep->d.err[i];
		} else {
			double mmodel = get_magnitude(a[index_mbar()+flt_id], a[index_beta()+flt_id], a[index_amp()], del_rho, del_tau);
			dy[i] = (cep->d.mag[i] - mmodel)/cep->d.err[i];
		}
	}

	return 0;
}

