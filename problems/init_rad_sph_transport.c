
#include "../decs.h"
#include "../constants.h"

void init_grid() {

	startx[0] = 0.;
	dx[0] = 5.e4;
	rtrans = 4.8e6;

	return;
}

void init_problem() {

	int g,ii,jj,kk,dd;
	double x,rho0,v0;

	gam = 5./3.;

	istep = 0;
 	t = 0.;

	rho0 = 1.e6;
	v0 = 0.;//0.05*CLIGHT;

	tmax = 4.e8/CLIGHT;

	t_next_restart = 1000.;
	t_next_log = 0.01;
	t_next_diag = 1000.;

	dt_restart = 1000.;
	dt_log = 0.01;
	dt_diag = 1000.;

	bound_flag[0][0] = REFLECT;
	bound_flag[0][1] = OUTFLOW;

	ZLOOP {
		x = r_of_x((ii+0.5)*dx[0] + startx[0]);
		NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho0;
		NDP_ELEM(sim.p,ii,jj,kk,UU) = 1.e6;//10*0.5*rho0*v0*v0;
		SLOOP {
			NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = 0.;
		}
		NDP_ELEM(sim.p,ii,jj,kk,U1) = v0;
		for(g=0;g<ngroups;g++) {
			NDP_ELEM(sim.p,ii,jj,kk,irad1+g) = exp(-pow(x-4.e8,2.)/(2.e7*2.e7))*100.;//1.e6*(-pow(g-nr1/2,2.)+pow(nr1/2+1.e-10,2.));  //exp(-pow(g-nr1/2,2.)/(10.*10.));
			DLOOP {
				NDP_ELEM(sim.p,ii,jj,kk,ifrad1+NDIM*g+dd) = -NDP_ELEM(sim.p,ii,jj,kk,irad1+g)*CLIGHT/ND_ELEM(geom,ii,jj,kk).scale[0][0];
			}
		}

		#if(NEUTRINO)
		NDP_ELEM(sim.p,ii,jj,kk,YE) = 0.5;
		#endif

	}


	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	return;
}
