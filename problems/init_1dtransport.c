
#include "../decs.h"
#include "../constants.h"

void init_grid() {

	startx[0] = 0.;
	dx[0] = CLIGHT/n1;

	return;
}

void init_problem() {

	int g,ii,jj,kk,dd;
	double T0, n0;

	gam = 5./3.;

	istep = 0;
 	t = 0.;

	dt_restart = 1000.;
	dt_log = 0.01;
	dt_diag = 1000.;

	bound_flag[0][0] = PROB;
	bound_flag[0][1] = OUTFLOW;

	ZLOOP {

		NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.;
		NDP_ELEM(sim.p,ii,jj,kk,UU) = 1.e6;
		DLOOP {
			NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = 0.;
		}
		for(g=0;g<ngroups;g++) {
			NDP_ELEM(sim.p,ii,jj,kk,irad1+g) = 0.;
			DLOOP {
				NDP_ELEM(sim.p,ii,jj,kk,ifrad1+NDIM*g+dd) = 0.;
			}
		}
	}

	reset_boundaries(sim.p);
	update_eos(sim.p, sim.eos);

	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	int g,dd;

	p[RHO] = 1.;
	p[UU] = 1.;
	DLOOP {
		p[U1+dd] = 0.;
	}

	for(g=0;g<ngroups;g++) {
		p[irad1+g] = 0.;
		DLOOP {
			p[ifrad1+NDIM*g+dd] = 0.;
		}
	}

	p[irad1] = 1.;
	p[ifrad1] = CLIGHT;

	return;
}
