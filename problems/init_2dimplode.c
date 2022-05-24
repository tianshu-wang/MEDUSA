
#include "../decs.h"

void init_grid() {

	startx[0] = 0.;
	startx[1] = 0.;
	dx[0] = 0.3/n1;
	dx[1] = 0.3/n2;

    periodic[0] = periodic[1] = periodic[2] = 0;

	return;
}

void init_problem() {

    int ii,jj,kk,vv,dd;
	double x,y;
	double pert_amp;

	gam = 1.4;

    VLOOP {
    DLOOP {
        bc[vv].lo[dd] = REFLECT;
        bc[vv].hi[dd] = REFLECT;
    }
    }

	ZLOOP {

		x = (ii+0.5)*dx[0] + startx[0];
		y = (jj+0.5)*dx[1] + startx[1];

		if(x+y > 0.15) {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.;
			NDP_ELEM(sim.p,ii,jj,kk,UU) = 1./(gam-1.);
		} else {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = 0.125;
			NDP_ELEM(sim.p,ii,jj,kk,UU) = 0.14/(gam-1.);
		}

		NDP_ELEM(sim.p,ii,jj,kk,U1) = 0.;
		NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.;
		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;

		//if(fabs(y-0.25) < dx[1] || fabs(y-0.75) < dx[1]) {
		//	NDP_ELEM(sim.p,ii,jj,kk,U2) = pert_amp*(2.*rand()/((double)RAND_MAX) - 1.);
		//}

	}

	return;
}

void prob_bounds(int i, int j, int k, double *p) {
}
