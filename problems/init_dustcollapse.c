#include "../decs.h"
#include "../constants.h"

static double rho0 = 1.0e9;
static double r0   = 6.5e8;
static double t0   = 1.0;
static double Mach = 2.0;
static double Mass0 = 0.0;

#if (ENFORCE_FLOORS!=TRUE)
#error This problem requires ENFORCE_FLOORS=TRUE
#endif

#if (EOS!=GAMMA_LAW)
#error This problem requires EOS=GAMMA_LAW
#endif

#if (GEOM!=SPHERICAL)
#error This problem requires GEOM=SPHERICAL
#endif

#if (GRAV!=SPHERICAL_MONOPOLE_GRAV)
#error This problem requires GRAV=SPHERICAL_MONOPOLE_GRAV
#endif

void init_grid()
{
	startx[0] = 0.0;
	dx[0]     = 5.0e4;
  #if (NDIM>1)
	startx[1] = -1.0;
	dx[1]     = 2.0/n2;
  #endif
  #if (NDIM==3)
  startx[2] = 0.0;
  dx[2]     = 2.0*M_PI/n3;
  #endif

  rtrans = rtrans_solve(dx[0], 7.0e8);  // approx transition radius from const to log spacing
  if (mpi_io_proc()) printf("[init_grid]:  rtrans=%e\n",rtrans);

	periodic[0] = 0;
	periodic[1] = 0;
	periodic[2] = 1;

	return;
}

void init_problem()
{
	double x[SPACEDIM];
  double r,v0,p0,rm,rp;
  int ii,jj,kk,vv;

  gam       = 5.0/3.0;
  v0        = r0/t0;
  p0        = rho0*SQR(v0/Mach)/gam;
	rho_floor = 1.0e-5*rho0;
  e_floor   = p0/(gam-1.0);

  VLOOP {
    bc[vv].lo[0] = SPHERICAL_ORIGIN;
    bc[vv].hi[0] = OUTFLOW;
    #if (NDIM>1)
    bc[vv].lo[1] = SPHERICAL_ORIGIN;
    bc[vv].hi[1] = SPHERICAL_ORIGIN;
    #endif
    #if (NDIM==3)
    bc[vv].lo[2] = PERIODIC;
    bc[vv].hi[2] = PERIODIC;
    #endif
  }
  bc[U1].hi[0] = PROB;

	ZLOOP {
    rm = r_of_x(startx[0] + (ii  )*dx[0]);
    rp = r_of_x(startx[0] + (ii+1)*dx[0]);
    if (rp < r0) {
      NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho0;
    } else if (rm < r0) {
      NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho0*(pow(r0,3) - pow(rm,3))/(pow(rp,3) - pow(rm,3));
    } else {
      NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho_floor;
    }    
    NDP_ELEM(sim.p,ii,jj,kk,UU ) = p0/(gam-1.0);
		NDP_ELEM(sim.p,ii,jj,kk,U1 ) = 0.0;
		NDP_ELEM(sim.p,ii,jj,kk,U2 ) = 0.0;
		NDP_ELEM(sim.p,ii,jj,kk,U3 ) = 0.0;
	}

	reset_boundaries(sim.p);
	update_eos(sim.p, sim.eos);

  fprintf(stderr,"myrank=%d, cs[0] = %e\n", myrank,NDP_ELEM(sim.eos,istart[0],istart[1],istart[2],CS));

	return;
}

void prob_bounds(int i, int j, int k, double *p)
{
  double r,vr;
  
  r     = r_of_x(beta0[i]/Gamma0[i]);
  vr    = -1.0/sqrt(r);
  p[U1] = vr/ND_ELEM(geom,i,j,k).scale[0][0];

	return;
}

void analysis_preloop()
{
  int ii,jj,kk;
  
  Mass0 = 0.0;
  ZLOOP { Mass0 += NDP_ELEM(sim.p,ii,jj,kk,RHO)*ND_ELEM(geom,ii,jj,kk).volume; }
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE, &Mass0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (myrank==0) printf("Mass=%e\n",Mass0);
  #endif

  return;
}

void analysis_inloop()
{
  // int ii,jj,kk;
  // double Mass = 0.0;
  //
  // ZLOOP { Mass += NDP_ELEM(sim.p,ii,jj,kk,RHO)*ND_ELEM(geom,ii,jj,kk).volume; }
  // #if (USE_MPI==TRUE)
  // MPI_Allreduce(MPI_IN_PLACE, &Mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // if (myrank==0) printf("Mass=%e, Relerr=%e\n",Mass,fabs(1.0-Mass/Mass0));
  // #endif

  return;
}

void analysis_postloop()
{
  int ii,jj,kk;
  double Mass = 0.0;
  
  ZLOOP { Mass += NDP_ELEM(sim.p,ii,jj,kk,RHO)*ND_ELEM(geom,ii,jj,kk).volume; }
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE, &Mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (myrank==0) printf("Mass=%e, Relerr=%e\n",Mass,fabs(1.0-Mass/Mass0));
  #endif

  return;
}