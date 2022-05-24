
#include "../../decs.h"
#include "../../constants.h"

void eos_p_cs(double *p, double *press, double *cs, double *Gamma) {

	*press = Kpoly*pow(p[RHO],gam);
	*cs = sqrt(gam*(*press)/p[RHO]);
	*Gamma = gam;

	return;
}

void eos_p(double *p, double *press) {

	*press = Kpoly*pow(p[RHO],gam);

	return;
}

double eos_cs(double *p) {

	double press,cs,Gamma;

	eos_p_cs(p, &press, &cs, &Gamma);

	return cs;
}

void eos_fill(double *p, double *eos, double avg_factor) {

	double press,cs,tp,Gamma;

	eos_p_cs(p, &press, &cs, &Gamma);
	tp = press/(p[RHO]*KBOLTZ)*MP;

	eos[PRESS] = press;
	eos[CS] = cs;
	eos[TEMP] = tp;
	eos[GAMMA] = Gamma;

}

void update_eos(double NDP_PTR p, double NDP_PTR eos, double avg_factor) {

	double Gamma;
	double press,cs;

	for(ii=-NG;ii<my_grid_dims[0]+NG;ii++) {
	#if(NDIM>1)
	for(jj=-NG;jj<my_grid_dims[1]+NG;jj++) {
	#endif
	#if(NDIM==3)
	for(kk=-NG;kk<my_grid_dims[2]+NG;kk++) {
	#endif
		eos_fill(ND_ELEM(p,ii,jj,kk), ND_ELEM(sim.eos,ii,jj,kk), avg_factor);
		//eos_p_cs(ND_ELEM(p,ii,jj,kk), &NDP_ELEM(eos,ii,jj,kk,PRESS), &NDP_ELEM(eos,ii,jj,kk,CS), &Gamma);
		/*eos_p_cs(ND_ELEM(p,ii,jj,kk), &press, &cs, &Gamma);
		if(avg_factor < 0.9) {
			NDP_ELEM(eos,ii,jj,kk,PRESS) += avg_factor*press;
			NDP_ELEM(eos,ii,jj,kk,CS) += avg_factor*cs;
		} else {
			NDP_ELEM(eos,ii,jj,kk,PRESS) = press;
			NDP_ELEM(eos,ii,jj,kk,CS) = cs;
		}*/
	#if(NDIM==3)
	}
	#endif
	#if(NDIM>1)
	}
	#endif
	}

	return;
}

