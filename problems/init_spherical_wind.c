#include "../decs.h"
#include "../constants.h"

static double drmin = 0.01;  // For Nr=128
static double rmax = 2.0;
static double alpha = 5.0;
static double a = 10.0;  // SET A
static double b = 0.0;  // SET A
// static double a = 16.0;  // SET B
// static double b = 0.5;  // SET B

static double rho_soln(double x);
static double vr_soln(double x);
static double p_soln(double x);
double rho_r2dr_dx(double x);
double vr_r2dr_dx(double x);
double p_r2dr_dx(double x);
double r2dr_dx(double x);

#if (ENFORCE_FLOORS!=FALSE)
#error This problem requires ENFORCE_FLOORS=FALSE
#endif

#if (GEOM!=SPHERICAL)
#error This problem requires GEOM=SPHERICAL
#endif

#if (GRAV!=NO_GRAV)
#error This problem requires GRAV=NO_GRAV
#endif

#if (EOS!=GAMMA_LAW)
#error This problem requires EOS=GAMMA_LAW
#endif


void init_grid()
{
	startx[0] = 0.0;
  dx[0]     = drmin*(128.0/n1);  // For RCOORD=SINH
  // dx[0]     = 2.0/n1;  // For RCOORD=UNIFORM
  #if (NDIM>1)
  startx[1] = -1.0;
  dx[1]     = 2.0/n2;
  #endif
  #if (NDIM==3)
  startx[2] = 0.0;
  dx[2]     = 2.0*M_PI/n3;
  #endif

  rtrans = rtrans_solve(dx[0], rmax);  // approx transition radius from const to log spacing
  if (myrank==0) printf("rtrans=%e\n",rtrans);

	periodic[0] = 0;
	periodic[1] = 0;
	periodic[2] = 1;

	return;
}

void init_problem()
{
	double x0lo,x0hi,r2dr_int,rho,vr,p;
  int ii,jj,kk,vv;

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

  gam = 5.0/3.0;

	ZLOOP {
    x0lo = startx[0] +  ii   *dx[0];
    x0hi = startx[0] + (ii+1)*dx[0];
    // if (jj==0) printf("r=%e, dr=%e\n",r_of_x(beta0[ii]/Gamma0[ii]),r_of_x(x0hi)-r_of_x(x0lo));
    
		r2dr_int = romb(x0lo, x0hi, r2dr_dx);
    rho = romb(x0lo, x0hi, rho_r2dr_dx)/r2dr_int;
    vr  = romb(x0lo, x0hi, vr_r2dr_dx)/r2dr_int;
    p   = romb(x0lo, x0hi, p_r2dr_dx)/r2dr_int;
    // if (jj==0) printf("ii=%d, rho=%e, vr=%e, p=%e, u=%e\n",ii,rho,vr,p,p/(gam-1.0));

    NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho;
    NDP_ELEM(sim.p,ii,jj,kk,UU ) = p/(gam-1.0);
		NDP_ELEM(sim.p,ii,jj,kk,U1 ) = vr/ND_ELEM(geom,ii,jj,kk).scale[0][0];
    NDP_ELEM(sim.p,ii,jj,kk,U2 ) = 0.0;
    NDP_ELEM(sim.p,ii,jj,kk,U3 ) = 0.0;
	}

	reset_boundaries(sim.p);
	update_eos(sim.p, sim.eos);

  fprintf(stderr,"myrank=%d, cs[0] = %e\n", myrank,NDP_ELEM(sim.eos,istart[0],istart[1],istart[2],CS));

	return;
}

static double rho_soln(double x)
{
  double tmp = 1.0/(1.0+alpha*t);
  double r = tmp*r_of_x(x);

  return pow(tmp,3)*(1.0 + exp(-pow(a,2)*pow(r-b,2)));
}

static double vr_soln(double x)
{
  return (alpha/(1.0+alpha*t))*r_of_x(x);
}

static double p_soln(double x)
{
  return pow(1.0/(1.0+alpha*t),3*gam)/gam;
}

double rho_r2dr_dx(double x)
{
  return rho_soln(x)*r2dr_dx(x);
}

double vr_r2dr_dx(double x)
{
  return vr_soln(x)*r2dr_dx(x);
}

double p_r2dr_dx(double x)
{
  return p_soln(x)*r2dr_dx(x);
}

void prob_bounds(int i, int j, int k, double *p)
{
	return;
}

void analysis_preloop()
{
  int ii,jj,kk;
  double x0lo,x0hi,r2dr_int,rho,rho_ref;
  double L1Err = 0.0;
  double vol_sum = 0.0;
  
  ZLOOP { vol_sum += ND_ELEM(geom,ii,jj,kk).volume; }
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE, &vol_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  
  ZLOOP {
    x0lo = startx[0] +  ii   *dx[0];
    x0hi = startx[0] + (ii+1)*dx[0];
  	r2dr_int = romb(x0lo, x0hi, r2dr_dx);
  
    rho = NDP_ELEM(sim.p,ii,jj,kk,RHO);
    rho_ref = romb(x0lo, x0hi, rho_r2dr_dx)/r2dr_int;
    
    L1Err += fabs(rho-rho_ref)*ND_ELEM(geom,ii,jj,kk).volume/vol_sum;
  }
  
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE, &L1Err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  
  if (myrank==0) printf("t=%e, n1=%d, L1Err=%e\n",t,n1,L1Err);
  
  return;
}

void analysis_inloop()
{
#if 0
  int ii,jj,kk;
  double x0lo,x0hi,r2dr_int,rho,rho_ref;
  double L1Err = 0.0;
  double vol_sum = 0.0;
  
  ZLOOP { vol_sum += ND_ELEM(geom,ii,jj,kk).volume; }
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE, &vol_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  
  ZLOOP {
    x0lo = startx[0] +  ii   *dx[0];
    x0hi = startx[0] + (ii+1)*dx[0];
  	r2dr_int = romb(x0lo, x0hi, r2dr_dx);
  
    rho = NDP_ELEM(sim.p,ii,jj,kk,RHO);
    rho_ref = romb(x0lo, x0hi, rho_r2dr_dx)/r2dr_int;
    
    L1Err += fabs(rho-rho_ref)*ND_ELEM(geom,ii,jj,kk).volume/vol_sum;
  }
  
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE, &L1Err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  
  if (myrank==0) printf("t=%e, n1=%d, L1Err=%e\n",t,n1,L1Err);
#endif
  
  return;
}

void analysis_postloop()
{
  int ii,jj,kk;
  double x0lo,x0hi,r2dr_int,rho,rho_ref;
  double L1Err = 0.0;
  double vol_sum = 0.0;
  
  ZLOOP { vol_sum += ND_ELEM(geom,ii,jj,kk).volume; }
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE, &vol_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  
  double L1Err_max = 0.0;
  double L1Err_global_max = 0.0;
  int max_loc = -1;
  ZLOOP {
    x0lo = startx[0] +  ii   *dx[0];
    x0hi = startx[0] + (ii+1)*dx[0];
  	r2dr_int = romb(x0lo, x0hi, r2dr_dx);
  
    rho = NDP_ELEM(sim.p,ii,jj,kk,RHO);
    rho_ref = romb(x0lo, x0hi, rho_r2dr_dx)/r2dr_int;
    
    L1Err += fabs(rho-rho_ref)*ND_ELEM(geom,ii,jj,kk).volume/vol_sum;
    if (fabs(rho-rho_ref)*ND_ELEM(geom,ii,jj,kk).volume/vol_sum > L1Err_max) {
      L1Err_max = fabs(rho-rho_ref)*ND_ELEM(geom,ii,jj,kk).volume/vol_sum;
      max_loc = ii;
    }
  }
  
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE, &L1Err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif  
  if (myrank==0) printf("t=%e, n1=%d, L1Err=%e\n",t,n1,L1Err);
  
  #if (USE_MPI==TRUE)
  MPI_Allreduce(&L1Err_max, &L1Err_global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  #endif
  if (L1Err_max==L1Err_global_max) printf("max_loc=%d\n",max_loc);
  
  
  return;  
}