#include "../decs.h"

void init_grid() {

	startx[0] = -1.2;
	startx[1] = -1.2;
	startx[2] = -1.2;
	dx[0] = 2.4/n1;
	dx[1] = 2.4/n2;
	dx[2] = 2.4/n3;

	periodic[0] = 0;
	periodic[1] = 0;
	periodic[2] = 0;	

	return;
}

void init_problem() {

	int ii,jj,kk,vv;
	double xl,xr;
	double rho_amb,e_amb,u_amb;
	double E_expl,R_expl;
	double x0[NDIM],x[NDIM];

	gam = 1.4;

	istep = 0;
 	t = 0.;
	dtmin = 1.e-15;

	rho_floor = 1.e-10;
	e_floor = 1.e-15;

	dt_restart = t_next_restart = 1000.;
	dt_log = t_next_log =  0.0005;
	dt_diag = t_next_diag = 0.0025;

	rho_amb = 1.;
	//e_amb = 1.6e-6/(gam-1.);
	//e_amb = 1.e-5/(gam-1.);
	//e_amb = 1./(gam-1.);
	e_amb = 1.e-12;
	u_amb = 0.;

	//E_expl = 0.244816;
	E_expl = 0.851072;
	R_expl = 0.075;
	//E_expl = 4./3.*M_PI*pow(4.e-3,3)*1.6e6/(gam-1.);

	VLOOP {
		bc[vv].lo[0] = OUTFLOW;
		#if(NDIM>1)
		bc[vv].lo[1] = OUTFLOW;
		bc[vv].hi[1] = OUTFLOW;
		#endif
		#if(NDIM==3)
		bc[vv].lo[2] = OUTFLOW;
		bc[vv].hi[2] = OUTFLOW;
		#endif
	}
	HLOOP {
		bc[vv].hi[0] = OUTFLOW;
	}

	srand(time(NULL));
	ZLOOP {
		NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho_amb;//*(ii+ 0.5)*dx[0];// + 0.5*sin(2.*M_PI*ii*dx[0]/0.3);
		NDP_ELEM(sim.p,ii,jj,kk,UU) = e_amb;
		NDP_ELEM(sim.p,ii,jj,kk,U1) = u_amb;
		NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.;//1.e-2*(2.*rand()/((double)RAND_MAX) - 1.);
		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;//10.;

	}

	int imax = R_expl/dx[0] + 1;
	double Edens = E_expl/(4./3.*M_PI*pow(R_expl,3.));
	double V0 = dx[0]*dx[1]*dx[2];

	double Etot = 0.;
	for(ii=istart[0];ii<istop[0];ii++) {
	for(jj=istart[1];jj<istop[1];jj++) {
	for(kk=istart[2];kk<istop[2];kk++) {
		// bottom left corner
		x0[0] = ii*dx[0]+startx[0];
		x0[1] = jj*dx[1]+startx[1];
		x0[2] = kk*dx[2]+startx[2];

		NDP_ELEM(sim.p,ii,jj,kk,UU) = 0.;
		NDP_ELEM(sim.p,ii,jj,kk,U1) = 0.;

		// subsample
		int nsamp = 8;
		double V=0.;
		double dV = V0/(nsamp*nsamp*nsamp);
		double E=0;
		for(int i=0;i<nsamp;i++) {
		for(int j=0;j<nsamp;j++) {
		for(int k=0;k<nsamp;k++) {
			x[0] = x0[0] + (i+0.5)*dx[0]/nsamp;
			x[1] = x0[1] + (j+0.5)*dx[1]/nsamp;
			x[2] = x0[2] + (k+0.5)*dx[2]/nsamp;
			double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
			if(r < R_expl) {
				E += dV*Edens;
				Etot += dV*Edens;
				V += dV;
			}
		}
		}
		}

		NDP_ELEM(sim.p,ii,jj,kk,UU) = E;
		NDP_ELEM(sim.p,ii,jj,kk,U1) = V;
		//if(V > V0) fprintf(stderr,"V > V0, wtf  %g %g\n", V, V0);

	}
	}
	}

	mpi_global_reduce(&Etot, 1);

	double Ecorr = E_expl/Etot;

	for(ii=istart[0];ii<istop[0];ii++) {
	for(jj=istart[1];jj<istop[1];jj++) {
	for(kk=istart[2];kk<istop[2];kk++) {

		NDP_ELEM(sim.p,ii,jj,kk,UU) = (NDP_ELEM(sim.p,ii,jj,kk,UU)*Ecorr + (V0-NDP_ELEM(sim.p,ii,jj,kk,U1))*e_amb)/V0;
		NDP_ELEM(sim.p,ii,jj,kk,U1) = 0.;
	}
	}
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
