
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

	rho_floor = 1.e-3;
	e_floor = 1.e-3;

	istep = 0;
 	t = 0.;

	t_next_dump = 0.005;
	t_next_restart = 1000.;
	t_next_log = 0.001;
	t_next_diag = 1000.;

	dt_dump = 0.005;
	dt_restart = 1000.;
	dt_log = 0.001;
	dt_diag = 1000.;

	tmax = 0.035;

	rhol = 5.99924;
	el = 460.894/(gam-1.);
	ul = 19.5975;

	rhor = 5.99242;
	er = 46.0950/(gam-1.);
	ur = -6.19633;

	fprintf(stderr,"my_grid_dims = %d\n", my_grid_dims[0]);

//	ZLOOP {
//		fprintf(stderr,"%d\n", ii);
//	}

	ZLOOP {
//		fprintf(stderr,"%d of %d\n", ii, my_grid_dims[0]);

		xl = ii*dx[0];
		xr = xl+dx[0];

		if(xr < 0.5) {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhol;
			NDP_ELEM(sim.p,ii,jj,kk,UU) = el;
			DLOOP {
				NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = ul;
			}
		} else if(xl < 0.5) {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhol*(0.5-xl)/dx[0] + rhor*(xr-0.5)/dx[0];
			NDP_ELEM(sim.p,ii,jj,kk,UU) = el*(0.5-xl)/dx[0] + er*(xr-0.5)/dx[0];
			DLOOP {
				NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = ul*(0.5-xl)/dx[0] + ur*(xr-0.5)/dx[0];
			}
		} else {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhor;
			NDP_ELEM(sim.p,ii,jj,kk,UU) = er;
			DLOOP {
				NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = ur;
			}
		}

	}

	reset_boundaries(sim.p);
	update_eos(sim.p, sim.eos);

	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	return;
}
