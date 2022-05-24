#include "../decs.h"
#include "../constants.h"

#if (GRAV!=USER_GRAV)
#error This problem must be compiled with GRAV=USER_GRAV
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

#if (NDIM!=3)
#error This problem must be compiled with NDIM=3
#endif

const double rmax = 1.0;  // max radius of sphere
const double Mach_max = 1.0;  // max Mach number of flow
// const double Mach_max = 10.0;  // max Mach number of flow
const double alpha = 0.0;  // alpha = {0, pi/2} for {z, x}-axis rotation
// const double alpha = 0.05*M_PI;  // alpha = {0, pi/2} for {z, x}-axis rotation
// const double alpha = 0.5*M_PI;  // alpha = {0, pi/2} for {z, x}-axis rotation
const double Omega0 = 2.0*M_PI;  // rotational velocity
const double rho0 = 1.0;  // density maximum
static double u0;  // constant internal energy
static double p0;
  
double vtheta(double r, double theta, double phi);
double vphi(double r, double theta, double phi);
double sth(double x);
double cthdth_dx(double x);
double sthdth_dx(double x);
double v2(int ii, int jj, int kk);
double v3(int ii, int jj, int kk);
void calculate_total_mass();

void init_grid()
{
	startx[0] = 0.0;
  startx[1] = -1.0;
  startx[2] = 0.0;
  dx[0]     = 0.01;
	dx[1]     = 2.0/n2;
  dx[2]     = 2.0*M_PI/n3;

  rtrans = rtrans_solve(dx[0], rmax);
  if (mpi_io_proc()) fprintf(stderr,"rtrans = %g\n", rtrans);

	periodic[0] = 0;
	periodic[1] = 0;
	periodic[2] = 1;

	return;
}

void init_problem()
{
	double r,th,ph,cs;
  int ii,jj,kk,vv,s,lev;

	VLOOP {
    bc[vv].lo[0] = SPHERICAL_ORIGIN;
    bc[vv].hi[0] = PROB;
    bc[vv].lo[1] = SPHERICAL_ORIGIN;
    bc[vv].hi[1] = SPHERICAL_ORIGIN;
		bc[vv].lo[2] = PERIODIC;
		bc[vv].hi[2] = PERIODIC;
	}

  gam = 5.0/3.0;
  cs = rmax*Omega0/Mach_max;
  rho_floor = 1.0e-3*rho0;
  p0 = rho0*SQR(cs)/gam;
  u0 = p0/(gam-1.0);
  e_floor = p0/(gam-1.0);

	ZLOOP {
    NDP_ELEM(sim.p,ii,jj,kk,RHO) = rho0;
    NDP_ELEM(sim.p,ii,jj,kk,UU ) = u0;
    NDP_ELEM(sim.p,ii,jj,kk,U1 ) = 0.0;
    NDP_ELEM(sim.p,ii,jj,kk,U2 ) = v2(ii,jj,kk);
    NDP_ELEM(sim.p,ii,jj,kk,U3 ) = v3(ii,jj,kk);
	}

	return;
}


void prob_bounds(int i, int j, int k, double *p)
{
  p[RHO] = rho0;
  p[UU ] = u0;
  p[U1 ] = 0.0;
  p[U2 ] = v2(i,j,k);
  p[U3 ] = v3(i,j,k);
  
	return;
}

double vtheta(double r, double theta, double phi)
{
  return r*Omega0*sin(phi)*sin(alpha);
}

double vphi(double r, double theta, double phi)
{
  return r*Omega0*(sin(theta)*cos(alpha) + cos(phi)*cos(theta)*sin(alpha));
}

double sth(double x)
{
	return sin(th_of_x(x));
}

double cthdth_dx(double x)
{
	return cos(th_of_x(x))*dth_dx(x);
}

/* Volume-averaged contravariant components of coordinate velocities v^2 and v^3
 * To obtain these, compute v_theta and v_phi, divide by sqrt(g_22) and sqrt(g_33), 
 * respectively, then take volume averages over sqrt(|g|) dx^1 dx^2 dx^3
 */
double v2(int ii, int jj, int kk)
{
  double x1lo,x1hi,vthf,x2lo,x2hi,vphf;

  x1lo = startx[1] +  jj   *DJS(ii)*dx[1];
  x1hi = startx[1] + (jj+1)*DJS(ii)*dx[1];
  vthf = romb(x1lo,x1hi,sthdth_dx);
  x2lo = startx[2] +  kk   *DKS(ii,jj)*dx[2];
  x2hi = startx[2] + (kk+1)*DKS(ii,jj)*dx[2];
  vphf = x2hi - x2lo;

  return Omega0*sin(alpha)*(cos(x2lo)-cos(x2hi))/(vthf*vphf)*romb(x1lo,x1hi,sth);
}

double v3(int ii, int jj, int kk)
{
  double x1lo,x1hi,vthf,x2lo,x2hi,vphf;

  x1lo = startx[1] +  jj   *DJS(ii)*dx[1];
  x1hi = startx[1] + (jj+1)*DJS(ii)*dx[1];
  vthf = romb(x1lo,x1hi,sthdth_dx);
  x2lo = startx[2] +  kk   *DKS(ii,jj)*dx[2];
  x2hi = startx[2] + (kk+1)*DKS(ii,jj)*dx[2];
  vphf = x2hi - x2lo;

  return Omega0*(cos(alpha) + sin(alpha)*(sin(x2hi)-sin(x2lo))/(vthf*vphf)*romb(x1lo,x1hi,cthdth_dx));
}

/* Compute the covariant components of the gravitational acceleration using the 
 * connection coefficients and stress tensor.  To obtain this, we require that 
 * the RHS of the equation written without covariant derivatives be 0; hence, 
 * \Gamma^l_{mn} T^m_l + \rho g_n = 0, where T^m_l = \rho v^m v_l + \delta^m_l P,
 * and \delta^m_l is the Kronecker symbol.
 */
void grav_accel(int i, int j, int k, double *g)
{
  int l,m,n;
  double vcon[3],vcov[3],T[3][3];
  
  vcon[0] = 0.0;
  vcon[1] = v2(i,j,k);
  vcon[2] = v3(i,j,k);
  geom_lower(vcon,vcov,ND_ELEM(geom,i,j,k).gcov[0]);
  
  for (l=0; l<3; l++) {
    for (m=0; m<3; m++) {
      T[m][l] = rho0*vcon[m]*vcov[l];
    }
    T[l][l] += p0;
  }
  
  memset(g, 0, 3*sizeof(double));  
  for (l=0; l<3; l++) {
    for (m=0; m<3; m++) {
      for (n=0; n<3; n++) {
        g[n] += -ND_ELEM(geom,i,j,k).conn[l][m][n]*T[m][l]/rho0;
      }
    }
  }

  // Raise the components, because a contravariant vector is expected!
  geom_raise(g,g,ND_ELEM(geom,i,j,k).gcon);
}

void analysis_preloop()
{
  calculate_total_mass();
}

void analysis_inloop()
{
  calculate_total_mass();
}

void analysis_postloop()
{
  calculate_total_mass();
}

// void calculate_total_mass()
// {
//   int ii,jj,kk;
//   double total_mass = 0.0;
//
//   ZLOOP { total_mass += NDP_ELEM(sim.p,ii,jj,kk,RHO)*ND_ELEM(geom,ii,jj,kk).volume; }
//   #if (USE_MPI==TRUE)
//   MPI_Allreduce(MPI_IN_PLACE,&total_mass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//   #endif
//   if (mpi_io_proc()) printf("Total mass = %1.15e\n",total_mass);
//
//   return;
// }