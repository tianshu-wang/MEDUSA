
#include "../decs.h"
#include "../constants.h"

void init_grid() {

	startx[0] = CLIGHT;
	//dx[0] = CLIGHT/n1;
	//startx[0] = 0.;
	dx[0] = CLIGHT/n1;
	rtrans = 100.*CLIGHT;//30.*CLIGHT/2.;

	startx[1] = -1.;
	dx[1] = 2./n2;

	return;
}

void init_problem() {

	int g,ii,jj,kk,dd,vv;
	double T0, n0;

	gam = 5./3.;

	istep = 0;
 	t = 0.;

	dt_restart = 1000.;
	dt_log = 0.01;
	dt_diag = 1000.;

	//bound_flag[0][0] = PROB;
	//bound_flag[0][1] = OUTFLOW;
	//bound_flag[0][0] = SPHERICAL_ORIGIN;
	//bound_flag[0][1] = PROB;

	VLOOP {
		bc[vv].lo[0] = PROB;
		bc[vv].hi[0] = OUTFLOW;
		#if(NDIM>1)
		bc[vv].lo[1] = REFLECT;
		bc[vv].hi[1] = REFLECT;
		#endif
	}


	ZLOOP {

		NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.;
		NDP_ELEM(sim.p,ii,jj,kk,UU) = 1.e6;
		SLOOP {
			NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = 0.;
		}
		#if(NEUTRINO==TRUE)
		NDP_ELEM(sim.p,ii,jj,kk,YE) = 0.5;
		#endif
		for(g=0;g<ngroups;g++) {
			NDP_ELEM(sim.p,ii,jj,kk,irad1+g) = 0.;
			DLOOP {
				NDP_ELEM(sim.p,ii,jj,kk,ifrad1+NDIM*g+dd) = 0.;
			}
		}
	}

	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	int g,dd;

	p[RHO] = 1.;
	p[UU] = 1.e6;
	DLOOP {
		p[U1+dd] = 0.;
	}
	#if(NEUTRINO==TRUE)
	p[YE] = 0.5;
	#endif

	for(g=0;g<ngroups;g++) {
		p[irad1+g] = 1.;// + 0.5*sin(2.*M_PI*g/(ngroups-1));
		double th = th_of_x((j+0.5)*dx[1] + startx[1]);

		// works
		p[ifrad1+NDIM*g] = chat*fabs(cos(th));
		p[ifrad1+NDIM*g+1] = chat*sin(th)*SIGN(th-M_PI/2.)/ND_ELEM(geom,i,j,k).scale[0][1];

		// works
		//p[ifrad1+NDIM*g] = chat*sin(th);
		//p[ifrad1+NDIM*g+1] = chat*cos(th)/ND_ELEM(geom,i,j,k).scale[0][1];

		// works
		//DLOOP {
		//	p[ifrad1+NDIM*g+dd] = 0.;//chat;//*(1 + 0.5*sin(2.*M_PI*g/(ngroups-1)));
		//}
		//p[ifrad1+NDIM*g] = chat;
	}

	//p[irad1] = 1.;//*CLIGHT*CLIGHT/(rcenter[i]*rcenter[i]);
	//p[ifrad1] = CLIGHT;// * CLIGHT*CLIGHT/(rcenter[i]*rcenter[i]);

	return;
}
