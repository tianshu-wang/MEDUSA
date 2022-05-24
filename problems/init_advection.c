#include "../decs.h"

void init_grid() {

	startx[0] = 0.;
	//startx[1] = 0.;
	startx[1] = 0.;
	startx[2] = 0.;
	//startx[2] = 0.;
	//dx[0] = 0.05/4.;
	//dx[0] = 0.705/n1;
	//dx[0] = 0.002;//0.175/n1;
    dx[0] = 1./n1;
	//dx[0] = 1.2/n1;
	//dx[1] = M_PI/n2;
	dx[1] = 1./n2;
	dx[2] = 1./n3;
	//rsparse_fact = 4.;
	//rsparse_trans = 0.075/2.;

	periodic[0] = 1;
	periodic[1] = 1;
	periodic[2] = 1;	

	return;
}

void init_problem() {

	int ii,jj,kk,dd,vv;
	double xl,xr;
	double rho_amb,e_amb,u_amb;

	gam = 1.4;

	istep = 0;
 	t = 0.;
	dtmin = 1.e-15;

	rho_floor = 1.e-10;
	e_floor = 1.e-15;

	rho_amb = 1.;
	//e_amb = 1.6e-6/(gam-1.);
	//e_amb = 1.e-5/(gam-1.);
	//e_amb = 1./(gam-1.);
	e_amb = 1e4*rho_amb/(gam*(gam-1));
	u_amb = 0.;

	VLOOP {
    DLOOP {
		bc[vv].lo[dd] = PERIODIC;
        bc[vv].hi[dd] = PERIODIC;
    }
	}

    if(restart_from_hdf) {
        restart_read();
        return;
    }

	//srand(time(NULL));
	ZLOOP {
        double x = startx[0] + (ii+0.5)*dx[0];
        if(x > 0.3 && x < 0.7) {
            NDP_ELEM(sim.p,ii,jj,kk,RHO) = 5*rho_amb;
        } else {
    		NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho_amb;//*(ii+0.5)*dx[0];// + 0.5*sin(2.*M_PI*ii*dx[0]/0.3);
        }
		NDP_ELEM(sim.p,ii,jj,kk,UU) = e_amb;
		NDP_ELEM(sim.p,ii,jj,kk,U1) = 1.;
		NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.;
		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;//10.;
	}
	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	int dd;
	p[RHO] = 1.;
	p[UU] = 1.e-12;
	SLOOP p[U1+dd] = 0.;

	return;
}
