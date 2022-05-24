#include "../decs.h"

void init_grid() {

	startx[0] = 0.;
	//startx[1] = 0.;
	startx[1] = -1.;
	startx[2] = 0.;
	//startx[2] = 0.;
	//dx[0] = 0.05/4.;
	//dx[0] = 0.705/n1;
	//dx[0] = 0.002;//0.175/n1;
    dx[0] = 0.4*1.2/n1;
	//dx[0] = 1.2/n1;
	//dx[1] = M_PI/n2;
	dx[1] = 2./n2;
	dx[2] = 2.*M_PI/n3;
	//rsparse_fact = 4.;
	//rsparse_trans = 0.075/2.;
	rtrans = rtrans_solve(dx[0], 1.2);
	if(mpi_io_proc()) fprintf(stderr,"rtrans = %g\n", rtrans);

	r_full_res = 0.02;

	//fprintf(stderr,"r comp   %g %g\n", r_of_x(dx[0]), rtrans*sinh(dx[0]/rtrans));

	periodic[0] = 0;
	periodic[1] = 0;
	periodic[2] = 1;	

	return;
}

void init_problem() {

	int ii,jj,kk,vv;
	double xl,xr;
	double rho_amb,e_amb,u_amb;
	double E_expl,R_expl;

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
	e_amb = 1.e-12;
	u_amb = 0.;

	//E_expl = 0.244816;
	E_expl = 0.851072;
	R_expl = 0.075;
	//E_expl = 4./3.*M_PI*pow(4.e-3,3)*1.6e6/(gam-1.);

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
	HLOOP {
		bc[vv].hi[0] = PROB;
	}

    if(restart_from_hdf) {
        restart_read();
        return;
    }

	//srand(time(NULL));
	ZLOOP {
		NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho_amb;//*(ii+0.5)*dx[0];// + 0.5*sin(2.*M_PI*ii*dx[0]/0.3);
		NDP_ELEM(sim.p,ii,jj,kk,UU) = e_amb;
		NDP_ELEM(sim.p,ii,jj,kk,U1) = u_amb;
		NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.;//1.e-4*(2*fornax_rand() - 1)/(r_of_x(startx[0]+(ii+0.5)*dx[0])*sin(th_of_x(startx[1]+(jj+0.5)*dx[1])));
		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;//10.;


/*
		double r = r_of_x((ii+0.5)*dx[0]);
		double th = (jj+0.5)*dx[1] + startx[1];
		double x = r*sin(th);
		double y = r*cos(th);
		double dist = sqrt(pow(x-0.075,2) + pow(y-0.0,2));
		double uring = 10.*E_expl/(2.*M_PI*M_PI*0.075*0.005*0.005);
		if(dist < 0.01) {
			NDP_ELEM(sim.p,ii,jj,kk,UU) = uring*exp(-pow(dist/0.01,2));
		}
		//NDP_ELEM(sim.p,ii,jj,kk,YE) = jj;
*/
	}

	double Edens = E_expl/(4./3.*M_PI*pow(R_expl,3.));

	if(redge[istart[0]] < R_expl) {
		for(ii=istart[0];ii<istop[0];ii++) {
			if(redge[ii+1] < R_expl) {
				#if(NDIM==1)
				NDP_ELEM(sim.p,ii,0,0,UU) =  Edens;
				#endif
				#if(NDIM>1)
				for(jj=istart[1];jj<istart[1] + my_grid_dims[1]/spider_fact[ii];jj++) {
				#if(NDIM==3)
				for(kk=istart[2];kk<istart[2] + my_grid_dims[2]/spider_fact[ii];kk++) {
				#endif
					NDP_ELEM(sim.p,ii,jj,kk,UU) = Edens;
				#if(NDIM==3)
				}
				#endif
				}
				#endif
			} else if(redge[ii] < R_expl) {
				Edens *= (pow(R_expl,3.)-pow(redge[ii],3.))/(pow(redge[ii+1],3.)-pow(redge[ii],3.));
				#if(NDIM==1)
				NDP_ELEM(sim.p,ii,0,0,UU) =  Edens;
				#endif
				#if(NDIM>1)
				for(jj=istart[1];jj<istart[1]+my_grid_dims[1]/spider_fact[ii];jj++) {
				#if(NDIM==3)
				for(kk=istart[2];kk<istart[2]+my_grid_dims[2]/spider_fact[ii];kk++) {
				#endif
					NDP_ELEM(sim.p,ii,jj,kk,UU) = Edens;
				#if(NDIM==3)
				}
				#endif
				}
				#endif
				break;
			} else {
				break;
			}				
		}
	}

/*	if(istart[0]<=300 && istop[0]>300) {
		if(istart[1] <= 128 && istop[1]>128) {
			NDP_ELEM(sim.p,300,128,0,UU) = E_expl/ND_ELEM(geom,300,128,0).volume;
		}
	}
*/
	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	int dd;
	p[RHO] = 1.;
	p[UU] = 1.e-12;
	SLOOP p[U1+dd] = 0.;

	return;
}
