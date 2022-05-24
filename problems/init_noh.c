
#include "../decs.h"

void init_grid() {

	startx[0] = 0.;
	startx[1] = -1.;
    startx[2] = 0.;

    
	dx[0] = 0.15/n1;
	dx[1] = 2.0/n2;
    dx[2] = 2*M_PI/n3;

    rtrans = rtrans_solve(dx[0], 1.0);
	if(mpi_io_proc()) fprintf(stderr,"rtrans = %g\n", rtrans);

	r_full_res = 0.1;

	periodic[0] = 0;
	periodic[1] = 0;
	periodic[2] = 1;

	return;
}

void init_problem() {

	int ii,jj,kk,vv;

	gam = 5./3.;
	rho_floor = 1.e-6;
	e_floor = 1.e-7/(gam-1.);

	VLOOP {
		bc[vv].lo[0] = SPHERICAL_ORIGIN;
		#if(NDIM>1)
		bc[vv].lo[1] = REFLECT;
		bc[vv].hi[1] = REFLECT;
		#endif
		#if(NDIM==3)
		bc[vv].lo[2] = PERIODIC;
		bc[vv].hi[2] = PERIODIC;
		#endif
	}
	VLOOP {
		bc[vv].hi[0] = PROB;
	}

	ZLOOP {

		NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.;
		NDP_ELEM(sim.p,ii,jj,kk,UU) = 1.e-6/(gam-1.);
		NDP_ELEM(sim.p,ii,jj,kk,U1) = -1./dr_dx((ii+0.5)*dx[0] + startx[0]);;
		NDP_ELEM(sim.p,ii,jj,kk,U2) = 1.e-8*(2*fornax_rand() - 1);
		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;
	}


	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	p[RHO] = pow(1. + t/r_of_x((i+0.5)*dx[0] + startx[0]),2.);
	p[U1] = -1./dr_dx((i+0.5)*dx[0] + startx[0]);
	p[U2] = 0.;
	p[U3] = 0.;
	p[UU] = 1.e-6/(gam-1.);

	return;
}
