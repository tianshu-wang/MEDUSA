
#include "../decs.h"

void init_grid() {

	startx[0] = 0.;
	dx[0] = 1./n1;

	return;
}

void init_problem() {

	int ii,jj,kk,dd,vv;
	double xl,xr;
	double rhol,el,ul;
	double rhor,er,ur;
	double rhoc,ec,uc;

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

	tmax = 0.038;

	rhol = 1.;
	el = 1000./(gam-1.);
	ul = 0.;

	rhoc = 1.;
	ec = 0.01/(gam-1.);
	uc = 0.;

	rhor = 1.;
	er = 100./(gam-1.);
	ur = 0.;

	// set the boundary conditions
	VLOOP {
		bc[vv].lo[0] = REFLECT;
        bc[vv].hi[0] = REFLECT;
    }

//	ZLOOP {
//		fprintf(stderr,"%d\n", ii);
//	}

	ZLOOP {
//		fprintf(stderr,"%d of %d\n", ii, my_grid_dims[0]);

		xl = ii*dx[0];
		xr = xl+dx[0];

		if(xr <= 0.1) {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhol;
			NDP_ELEM(sim.p,ii,jj,kk,UU) = el;
			NDP_ELEM(sim.p,ii,jj,kk,U1) = ul;
		}

		if(xl < 0.1 && xr > 0.1) {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhol*(0.1-xl)/dx[0] + rhoc*(xr-0.1)/dx[0];
			NDP_ELEM(sim.p,ii,jj,kk,UU) = el*(0.1-xl)/dx[0] + ec*(xr-0.1)/dx[0];
			NDP_ELEM(sim.p,ii,jj,kk,U1) = ul*(0.1-xl)/dx[0] + uc*(xr-0.1)/dx[0];
		}

		if(xl >= 0.1 && xr <= 0.9) {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhoc;
			NDP_ELEM(sim.p,ii,jj,kk,UU) = ec;
			NDP_ELEM(sim.p,ii,jj,kk,U1) = uc;
		}

		if(xl <= 0.9 && xr > 0.9) {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhoc*(0.9-xl)/dx[0] + rhor*(xr-0.9)/dx[0];
			NDP_ELEM(sim.p,ii,jj,kk,UU) = ec*(0.9-xl)/dx[0] + er*(xr-0.9)/dx[0];
			NDP_ELEM(sim.p,ii,jj,kk,U1) = uc*(0.9-xl)/dx[0] + ur*(xr-0.9)/dx[0];
		}

		if(xl >= 0.9) {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhor;
			NDP_ELEM(sim.p,ii,jj,kk,UU) = er;
			NDP_ELEM(sim.p,ii,jj,kk,U1) = ur;
		}
			
		NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.;
		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;	

	}

//	reset_boundaries(sim.p);
//	update_eos(sim.p, sim.eos);

	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	return;
}
