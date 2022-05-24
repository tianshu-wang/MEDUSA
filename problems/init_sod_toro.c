
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

	dt_dump = 0.001;
	dt_restart = 1000.;
	dt_log = 0.01;
	dt_diag = 1000.;

	tmax = 0.2;

	rhol = 1.;
	el = 1./(gam-1.);
	ul = 0.75;

	rhor = 0.125;
	er = 0.1/(gam-1.);
	ur = 0.;

	fprintf(stderr,"my_grid_dims = %d\n", my_grid_dims[0]);

//	ZLOOP {
//		fprintf(stderr,"%d\n", ii);
//	}

	ZLOOP {
//		fprintf(stderr,"%d of %d\n", ii, my_grid_dims[0]);

		xl = ii*dx[0];
		xr = xl+dx[0];

		if(xr < 0.3) {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhol;
			NDP_ELEM(sim.p,ii,jj,kk,UU) = el;
			NDP_ELEM(sim.p,ii,jj,kk,U1) = ul;
		} else if(xl < 0.3) {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhol*(0.3-xl)/dx[0] + rhor*(xr-0.3)/dx[0];
			NDP_ELEM(sim.p,ii,jj,kk,UU) = el*(0.3-xl)/dx[0] + er*(xr-0.3)/dx[0];
			NDP_ELEM(sim.p,ii,jj,kk,U1) = ul*(0.3-xl)/dx[0] + ur*(xr-0.3)/dx[0];
		} else {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhor;
			NDP_ELEM(sim.p,ii,jj,kk,UU) = er;
			NDP_ELEM(sim.p,ii,jj,kk,U1) = ur;
		}

	}

	reset_boundaries(sim.p);
	update_eos(sim.p, sim.eos);

	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	return;
}
