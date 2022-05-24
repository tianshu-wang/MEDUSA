
#include "../decs.h"

void init_grid() {

	startx[0] = 0.;
	dx[0] = 1./n1;

	return;
}

void init_problem() {

	int ii,jj,kk,dd;
	double xl,xr;
	double rhol,el,ul;
	double rhor,er,ur;

	gam = 1.4;

	istep = 0;
 	t = 0.;

	t_next_dump = 0.05;
	t_next_restart = 1000.;
	t_next_log = 0.01;
	t_next_diag = 1000.;

	dt_dump = 0.05;
	dt_restart = 1000.;
	dt_log = 0.01;
	dt_diag = 1000.;

	tmax = 0.25;

	ZLOOP {

		NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.;
		NDP_ELEM(sim.p,ii,jj,kk,UU) = 1. + 1.*(ii%2);
		DLOOP {
			NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = 0.;
		}
	}

	reset_boundaries(sim.p);
	update_eos(sim.p, sim.eos);

	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	return;
}
