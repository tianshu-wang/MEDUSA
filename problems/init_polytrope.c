#include "../decs.h"
#include "../constants.h"

void init_grid() {

	startx[0] = 0.;
	startx[1] = 0.;
	dx[0] = 5.e4;
	dx[1] = M_PI/N2;
	rtrans = 4.8e6;

	return;
}

void init_problem(char *name) {

	double xl,xr,xi,alpha;

	gam = 2.;
	Kpoly = 1.55e5;
	alpha = sqrt(2.*Kpoly/(4.*M_PI*GNEWT));

	istep = 0;
 	t = 0.;
	dtmin = 1.e-15;

	rho_floor = 0.;
	e_floor = 0.;

	dt_dump = t_next_dump = 1.e-2;
	dt_restart = t_next_restart = 1000.;
	dt_log = t_next_log =  0.0005;
	dt_diag = t_next_diag = 0.0025;

	tmax = 1.;

	fprintf(stdout,"dx[0] = %g\n", dx[0]);

	bound_flag[0][0] = REFLECT;
	bound_flag[0][1] = PROB;
	bound_flag[1][0] = REFLECT;
	bound_flag[1][1] = REFLECT;

	ZLOOP {
		xi = rcenter[ii]/alpha;
		if(xi < M_PI) {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = 2.e14*sin(xi)/xi + 1.e-5;	// + 1 so it doesn't ever go to zero
		} else {
			NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.e-5;// * (1910140.114312/rcenter[ii]);// + GNEWT*1.7679355942703e33/(Kpoly*gam)*(1./rcenter[ii] - 1./1910140.114312);
		}
		NDP_ELEM(sim.p,ii,jj,kk,UU) = NDP_ELEM(sim.p,ii,jj,kk,RHO);
		NDP_ELEM(sim.p,ii,jj,kk,U1) = 0.;
		NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.;
		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;
	}


	reset_boundaries(sim.p);
	update_eos(sim.p, sim.eos, 1.);
	
	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	p[RHO] = 1.e-5;
	p[UU] = 1.e-5;
	p[U1] = 0.;
	p[U2] = 0.;
	p[U3] = 0.;

	return;
}
