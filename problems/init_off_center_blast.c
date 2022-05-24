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

void init_grid()
{
	startx[0] = 0.0;
	startx[1] = -1.0;
	startx[2] = 0.0;
  dx[0]     = 0.4*1.2/n1;
  dx[1]     = 2.0/n2;
  dx[2]     = 2.0*M_PI/n3;
  
  rtrans = rtrans_solve(dx[0], 1.2);
  if (mpi_io_proc()) fprintf(stderr,"rtrans = %g\n", rtrans);

	periodic[0] = 0;
	periodic[1] = 0;
	periodic[2] = 1;

	return;
}

void init_problem()
{
  int ii,jj,kk,vv,iexp,jexp,kexp,jp,kp;
  double xexp,yexp,zexp,rexp,r,th,ph,xx,yy,zz;
	double xl,xr;
	double rho_amb,e_amb,u_amb,Edens;
	double E_expl,R_expl;
  double x[3];

	gam = 1.4;

	rho_amb = 1.0;
	e_amb   = 1.0e-12;
	u_amb   = 0.0;

	E_expl = 0.851072;
	R_expl = 0.075;
 
  r  = 0.12;
  // th = 0.25*M_PI;
  th = 0.5*M_PI;
  // r  = 0.5;
  // th = 0.1*M_PI;
  ph = 0.0;
  
  xexp = r*sin(th)*cos(ph);
  yexp = r*sin(th)*sin(ph);
  zexp = r*cos(th);
  if (myrank==0) printf("xexp=%e, yexp=%e, zexp=%e\n",xexp,yexp,zexp);
 
	VLOOP {
    bc[vv].lo[0] = SPHERICAL_ORIGIN;
    #if (NDIM>1)
    bc[vv].lo[1] = SPHERICAL_ORIGIN;
    bc[vv].hi[1] = SPHERICAL_ORIGIN;
    #endif
		#if (NDIM==3)

		bc[vv].lo[2] = PERIODIC;
		bc[vv].hi[2] = PERIODIC;
		#endif
	}
	HLOOP {
		bc[vv].hi[0] = PROB;
	}

  if (restart_from_hdf) {
    restart_read();
    return;
  }

	Edens = E_expl/(4.0/3.0*M_PI*pow(R_expl,3.0));

	ZLOOP {
		NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho_amb;
		NDP_ELEM(sim.p,ii,jj,kk,UU ) = e_amb;
		NDP_ELEM(sim.p,ii,jj,kk,U1 ) = u_amb;
		NDP_ELEM(sim.p,ii,jj,kk,U2 ) = 0.0;
		NDP_ELEM(sim.p,ii,jj,kk,U3 ) = 0.0;

    ijk_to_x(ii,jj,kk,x);
    r  = r_of_x(x[0]);
    th = th_of_x(x[1]);
    ph = x[2];
    xx = r*sin(th)*cos(ph);
    yy = r*sin(th)*sin(ph);
    zz = r*cos(th);
    
    if (SQR(xx-xexp) + SQR(yy-yexp) + SQR(zz-zexp) < SQR(R_expl)) {
      NDP_ELEM(sim.p,ii,jj,kk,UU) = Edens;
    }
  }

	return;
}

void prob_bounds(int i, int j, int k, double *p)
{
	int dd;
	
	p[RHO] = 1.0;
	p[UU ] = 1.0e-12;
	SLOOP { p[U1+dd] = 0.0; }

	return;
}

void analysis_preloop()
{
  return;
}

void analysis_inloop()
{
  // int ii,jj,kk,vv;
  //
  // // check for Inf and NaN
  // ZLOOP {
  //   VLOOP {
  //     if (vv==UU && ii<10 && isnan(NDP_ELEM(sim.p,ii,jj,kk,vv)))
  //       printf("[analysis_inloop]:  proc %d detected NaN in variable %d at zone (%d,%d,%d)\n",myrank,vv,ii,jj,kk);
  //     if (vv==UU && ii<10 && isinf(NDP_ELEM(sim.p,ii,jj,kk,vv)))
  //       printf("[analysis_inloop]:  proc %d detected Inf in variable %d at zone (%d,%d,%d)\n",myrank,vv,ii,jj,kk);
  //   }
  // }

  return;
}

void analysis_postloop()
{
  return;
}