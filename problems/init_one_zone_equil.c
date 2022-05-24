
#include "../decs.h"
#include "../constants.h"

void init_grid() {

	startx[0] = 0.;
	dx[0] = 0.01*CLIGHT/n1;
	rtrans = 100.*CLIGHT;

	return;
}

void init_problem() {

	int g,ii,jj,kk,dd;
	double T0, n0;

	gam = 5./3.;

	istep = 0;
 	t = 0.;

	//t_next_dump = 0.05;
	t_next_restart = 100000.;
	t_next_log = 0.01;
	t_next_diag = 100000.;

	//dt_dump = 0.05;
	dt_restart = 100000.;
	dt_log = 0.01;
	dt_diag = 100000.;

	T0 = 1.e6;
	n0 = HPLANCK*sqrt(T0)/((gam-1.)*M_PI*(6.8e-38/(4.*M_PI)))*10.;

	ZLOOP {

		NDP_ELEM(sim.p,ii,jj,kk,RHO) = n0*(MP+ME);
		NDP_ELEM(sim.p,ii,jj,kk,UU) = 2.*n0*KBOLTZ*T0/(gam-1.);
		DLOOP {
			NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = 0.;
		}
		NDP_ELEM(sim.p,ii,jj,kk,U1) = 5.e8;
		for(g=0;g<ngroups;g++) {
			NDP_ELEM(sim.p,ii,jj,kk,irad1+NDIM*g) = 0.;
			DLOOP {
				NDP_ELEM(sim.p,ii,jj,kk,irad1+NDIM*g+dd) = 0.;
			}
		}
	}

	//reset_boundaries(sim.p);
	//update_eos(sim.p, sim.eos);

	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	return;
}
