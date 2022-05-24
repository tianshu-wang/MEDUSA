#include "../decs.h"

#if (GRAV!=NO_GRAV)
#error This problem must be compiled with GRAV=NO_GRAV
#endif

#if (GEOM!=SPHERICAL)
#error This problem must be compiled with GEOM=SPHERICAL
#endif

#if (EOS!=GAMMA_LAW)
#error This problem must be compiled with EOS=GAMMA_LAW
#endif

#if (USE_EXT_SRC!=FALSE)
#error This problem must be compiled with USE_EXT_SRC=FALSE
#endif

#if (USE_AUX_SRC!=FALSE)
#error This problem must be compiled with USE_AUX_SRC=FALSE
#endif

void init_grid() {

	startx[0] = 0.;
	startx[1] = -1.;
	dx[0] = 0.4/n1;
	dx[1] = 2./n2;
	dx[2] = 2*M_PI/n3;	


	//dx[0] = 1.e4;
	rsparse_fact = 4.;
	rsparse_trans = 0.05;
	rtrans = rtrans_solve(dx[0], 0.8);

	periodic[0] = 0;
	periodic[1] = 0;
	periodic[2] = 1;

	return;
}

void init_problem() {

	int ii,jj,kk,vv;
	double xl,xr;
	double rho_amb,e_amb,u_amb;
	double E_expl;

	gam = 1.4;

	dtmin = 1.e-15;

	rho_floor = 1.e-10;
	e_floor = 1.e-10;

	dt_log = t_next_log =  0.0005;
	dt_diag = t_next_diag = 0.0025;


	rho_amb = 1.;
	e_amb = 1.e-5/(gam-1.);
	u_amb = 0.;

	E_expl = 1.;

	// set the boundary conditions
	VLOOP {
		bc[vv].lo[0] = REFLECT;
		#if(NDIM>1)
		bc[vv].lo[1] = REFLECT;
		bc[vv].hi[1] = REFLECT;
		#endif
		#if(NDIM==3)
		bc[vv].lo[2] = PERIODIC;
		bc[vv].hi[2] = PERIODIC;
		#endif
	}
	HLOOP {
		bc[vv].hi[0] = EXTRAP;
	}

	ZLOOP {
		NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho_amb;//*(ii+0.5)*dx[0];// + 0.5*sin(2.*M_PI*ii*dx[0]/0.3);
		NDP_ELEM(sim.p,ii,jj,kk,UU) = e_amb;
		NDP_ELEM(sim.p,ii,jj,kk,U1) = u_amb;
		NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.;
		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;
	}
	
	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	return;
}
