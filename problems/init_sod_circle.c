
#include "../decs.h"

void init_grid() {

	startx[0] = 0.;
	startx[1] = 0.;
    startx[2] = 0.;
	dx[0] = 2.0/n1;
	dx[1] = 2.0/n2;
    dx[2] = 1.0/n3;

    periodic[0] = periodic[1] = periodic[2] = 0;

	return;
}

void init_problem() {

    int ii,jj,kk,vv,dd;
	double x,y;
	double pert_amp;

	gam = 5./3.;

    VLOOP {
    DLOOP {
        bc[vv].lo[dd] = REFLECT;
        bc[vv].hi[dd] = REFLECT;
    }
    }

	ZLOOP {

		x = (ii+0.5)*dx[0] + startx[0];
		y = (jj+0.5)*dx[1] + startx[1];
        double r = sqrt(x*x + y*y);

        NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.;
        if(r < 1.) {
            NDP_ELEM(sim.p,ii,jj,kk,UU) = 1.e-6;
        } else {
            NDP_ELEM(sim.p,ii,jj,kk,UU) = 1.;
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
