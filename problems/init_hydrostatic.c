
#include "../decs.h"
#include "../constants.h"

#if (!(GRAV==SPHERICAL_MONOPOLE_GRAV || GRAV==SPHERICAL_MULTIPOLE_GRAV))
#error This problem must be compiled with GRAV=SPHERICAL_MONOPOLE_GRAV or GRAV=SPHERICAL_MULTIPOLE_GRAV
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

double ran2(long int *idum);

void init_grid()
{
	startx[0] = 0.0;
	startx[1] = -1.0;
  startx[2] = 0.0;
  dx[0]     = 5.0e4;
  dx[1]     = 2.0/n2;
  dx[2]     = 2.0*M_PI/n3;
  
  rtrans = rtrans_solve(dx[0], 1.0e9);
  if (mpi_io_proc()) fprintf(stderr,"rtrans = %g\n", rtrans);

	periodic[0] = 0;
	periodic[1] = 0;
	periodic[2] = 1;

	return;
}

static double rho0,alpha;

void init_problem()
{
  int vv,ii,jj,kk;
	double x[SPACEDIM];
  double* r_model=NULL;
  double* rho_model=NULL;
  double* uu_model=NULL;
  double* v_model=NULL;
  double r,rho,uu,v,dummy,del;
  char ch,dum[1024];
  int i,idum,lines;
  FILE *fp=NULL;

	gam       = 5.0/3.0;
	rho_floor = 1.0e-6;
	e_floor   = 1.0e-9/(gam-1.0);

	VLOOP {
    bc[vv].lo[0] = SPHERICAL_ORIGIN;
    #if(NDIM>1)
    bc[vv].lo[1] = SPHERICAL_ORIGIN;
    bc[vv].hi[1] = SPHERICAL_ORIGIN;
		#endif
		#if(NDIM==3)
		bc[vv].lo[2] = PERIODIC;
		bc[vv].hi[2] = PERIODIC;
		#endif
	}
	HLOOP {
    bc[vv].hi[0] = PROB;
	}

  alpha = -2.0;  // do not set alpha = -1.0
  rho0 = 1.0e14*pow(r_of_x(0.5*dx[0]),-alpha);
	if (mpi_io_proc()) fprintf(stderr,"rho0 = %g\n", rho0);

	ZLOOP {
    ijk_to_x(ii,jj,kk,x);
    r = r_of_x(x[0]);
		NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho0*pow(r,alpha);
    NDP_ELEM(sim.p,ii,jj,kk,UU ) = -2.0*M_PI*GNEWT*SQR(rho0)/((alpha+3.0)*(alpha+1.0))*pow(r,2.0*alpha+2.0)/(gam-1.0);
		NDP_ELEM(sim.p,ii,jj,kk,U1 ) = 0.0;
		NDP_ELEM(sim.p,ii,jj,kk,U2 ) = 0.0;
		NDP_ELEM(sim.p,ii,jj,kk,U3 ) = 0.0;
	}
  
	reset_boundaries(sim.p);
	update_eos(sim.p, sim.eos);

  long int rseed = -myrank;
  #if 0
  ZLOOP { NDP_ELEM(sim.p,ii,jj,kk,RHO) *= 1.0 + (perturb_level/100.0)*(2.0*ran2(&rseed) - 1.0); }
  #else
  ZLOOP {
    NDP_ELEM(sim.p,ii,jj,kk,U1) += NDP_ELEM(sim.eos,istart[0],istart[1],istart[2],CS)*(perturb_level/100.0)*(2.0*ran2(&rseed) - 1.0)/ND_ELEM(geom,ii,jj,kk).scale[0][0];
    NDP_ELEM(sim.p,ii,jj,kk,U2) += NDP_ELEM(sim.eos,istart[0],istart[1],istart[2],CS)*(perturb_level/100.0)*(2.0*ran2(&rseed) - 1.0)/ND_ELEM(geom,ii,jj,kk).scale[0][1];
    NDP_ELEM(sim.p,ii,jj,kk,U3) += NDP_ELEM(sim.eos,istart[0],istart[1],istart[2],CS)*(perturb_level/100.0)*(2.0*ran2(&rseed) - 1.0)/ND_ELEM(geom,ii,jj,kk).scale[0][2];
  }
  #endif

  fprintf(stderr,"myrank=%d, cs[0] = %e\n", myrank,NDP_ELEM(sim.eos,istart[0],istart[1],istart[2],CS));

	return;
}

void prob_bounds(int i, int j, int k, double *p)
{
	double r,x[SPACEDIM];

	ijk_to_x(i,j,k,x);
	r = r_of_x(x[0]);
	p[RHO] = rho0*pow(r,alpha);
  p[UU ] = -2.0*M_PI*GNEWT*SQR(rho0)/((alpha+3.0)*(alpha+1.0))*pow(r,2.0*alpha+2.0)/(gam-1.0);
	p[U1 ] = 0.0;
	p[U2 ] = 0.0;
	p[U3 ] = 0.0;

	return;
}