
#include "../decs.h"

void init_grid() {

	startx[0] = 0.0;
	startx[1] = 0.;
	dx[0] = 1.0/N1;
	dx[1] = M_PI/N2;

	return;
}

static double vy;

void init_problem() {

	double x,y,r,th,sth,cth;
	double vr,vth;
	double xb,yb,rb;

	gam = 5./3.;
	rho_floor = 1.e-6;
	e_floor = 1.e-9/(gam-1.);

	istep = 0;
 	t = 0.;

	dt_dump = 0.025;
	dt_restart = 1000.;
	dt_log = 0.01;
	dt_diag = 1000.;

	t_next_dump = dt_dump;
	t_next_restart = dt_restart;
	t_next_log = dt_log;
	t_next_diag = dt_diag;

	tmax = 1.5;

	bound_flag[0][0] = REFLECT;
	bound_flag[0][1] = PROB;
	bound_flag[1][0] = REFLECT;
	bound_flag[1][1] = REFLECT;

	xb = 0.3;
	yb = 0.75;

	vy = 0.;
	ZLOOP {
		r = (ii+0.5)*dx[0] + startx[0];
		th = (jj+0.5)*dx[1] + startx[1];
		sth = sin(th);
		cth = cos(th);
		x = r*sth;
		y = r*cth;
		vr = vy*y/r;
		vth = -vy*x/(r*r);
		rb = sqrt(pow(x-xb,2.) + pow(y-yb,2.));
		NDP_ELEM(sim.p,ii,jj,kk,RHO) = 2.-(atan((rb-0.1)/dx[0])+0.5*M_PI)/M_PI;
		//NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1. + 0.5*sin(5.*th)*r*r;
		//if(y < 0.25) NDP_ELEM(sim.p,ii,jj,kk,RHO) = 2.;
		//else NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.;
		NDP_ELEM(sim.p,ii,jj,kk,UU) = 10./(gam-1.);
		//if(x > 0.1) {
			NDP_ELEM(sim.p,ii,jj,kk,U1) = vr;
			NDP_ELEM(sim.p,ii,jj,kk,U2) = vth;
		//} else {
		//	NDP_ELEM(sim.p,ii,jj,kk,U1) = 0.;
		//	NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.;
		//}

		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;
	}

	reset_boundaries(sim.p);

	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	double x,y,r,th,sth,cth;

	p[RHO] = 1.;
	p[UU] = 10./(gam-1.);
	r = (i+0.5)*dx[0] + startx[0];
	th = (j+0.5)*dx[1] + startx[1];
	sth = sin(th);
	cth = cos(th);
	x = r*sth;
	y = r*cth;
	//if(x > 0.1) {
		p[U1] = vy*y/r;
		p[U2] = -vy*x/(r*r);
	//} else {
	//	p[U1] = 0.;
	//	p[U2] = 0.;
	//}
	p[U3] = 0.;

	return;
}
