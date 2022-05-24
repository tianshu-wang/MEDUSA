
#include "../decs.h"
#include "../constants.h"

void init_grid() {

	startx[0] = -1.;
	dx[0] = 2./n1;
	//rtrans = CLIGHT;

	return;
}

void init_problem() {

	int g,ii,jj,kk,dd,vv;
	double x,rho0,v0;

	gam = 5./3.;

	istep = 0;
 	t = 0.;

	rho0 = 1.e6;
	v0 = 0.;//5.e8;//0.;//0.05*CLIGHT;

	//tmax = 2.e-9;//1./v0;	// should advect back to where it started

	t_next_restart = 1000.;
	t_next_log = 0.01;
	t_next_diag = 1000.;

	dt_restart = 1000.;
	dt_log = 0.01;
	dt_diag = 1000.;

	// set the boundary conditions
	VLOOP {
		bc[vv].lo[0] = OUTFLOW;
		bc[vv].hi[0] = OUTFLOW;
	}

	/*for(vv=irad1;vv<irad1+(1+NDIM)*ngroups;vv++) {
		bc[vv].hi[0] = RADEXTRAP;
	}*/

	ZLOOP {
		x = (ii+0.5)*dx[0] + startx[0];
		NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho0;
		NDP_ELEM(sim.p,ii,jj,kk,UU) = 10*0.5*rho0*v0*v0 + rho0;
		SLOOP {
			NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = 0.;
		}
		NDP_ELEM(sim.p,ii,jj,kk,U1) = v0;
		for(g=0;g<ngroups;g++) {
			NDP_ELEM(sim.p,ii,jj,kk,irad1+g) = exp(-x*x/(0.1*0.1));//*(-pow(g-nr1/2,2.)+pow(nr1/2+1.e-10,2.));  //exp(-pow(g-nr1/2,2.)/(10.*10.));
			DLOOP {
				NDP_ELEM(sim.p,ii,jj,kk,ifrad1+NDIM*g+dd) = 2*CLIGHT/(3.*512.)*x/(0.1*0.1)*exp(-x*x/(0.1*0.1));
//0.;//NDP_ELEM(sim.p,ii,jj,kk,irad1+g)*NDP_ELEM(sim.p,ii,jj,kk,U1);
			}
		}

		#if(NEUTRINO)
		NDP_ELEM(sim.p,ii,jj,kk,YE) = 0.5;
		#endif

	}


	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	return;
}
