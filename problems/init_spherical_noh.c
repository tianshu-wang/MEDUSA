
#include "../decs.h"

void init_grid() {

	startx[0] = 0.;
	startx[1] = 0.;
	dx[0] = 1.0/n1;
	dx[1] = M_PI/n2;
	rtrans = 3.;

	return;
}

void init_problem() {

	int ii,jj,kk,vv;
	double x,y,r;
	double vr;

	gam = 5./3.;
	rho_floor = 1.e-6;
	e_floor = 1.e-7/(gam-1.);

	istep = 0;
 	t = 0.;

	dt_restart = 1000.;
	dt_log = 0.01;
	dt_diag = 1000.;

	t_next_dump = dt_dump;
	t_next_restart = dt_restart;
	t_next_log = dt_log;
	t_next_diag = dt_diag;


	// set the boundary conditions
	VLOOP {
		bc[vv].lo[0] = REFLECT;
		#if(NDIM>1)
		bc[vv].lo[1] = REFLECT;
		bc[vv].hi[1] = REFLECT;
		#endif
	}
	VLOOP {
		bc[vv].hi[0] = PROB;
	}

	vr = -1.;

	ZLOOP {

		NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.;
		NDP_ELEM(sim.p,ii,jj,kk,UU) = 1.e-6/(gam-1.);
		NDP_ELEM(sim.p,ii,jj,kk,U1) = vr/ND_ELEM(geom,ii,jj,kk).scale[0][0];
		NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.;
		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;
	}

	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	double x,y,r,rhat[SPACEDIM];

	p[RHO] = 1.;
	p[UU] = 1.e-6/(gam-1.);
	r = ijk_to_r(i,j,k,rhat);//(i+0.5)*dx[0] + startx[0];
	p[U1] = -1./rhat[0];
	//fprintf(stderr,"bound: %g %g\n", r, p[U1]);
	p[U2] = 0.;
	p[U3] = 0.;

	p[RHO] += t/r;
	p[RHO] *= p[RHO];

	return;
}
