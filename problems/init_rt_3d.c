#include "../decs.h"
#include "../constants.h"

/*
* Sets up a 3D KH test as in Phys Fluids 16:1668 (2004)
*/

#if (DO_HYDRO!=TRUE)
#error This problem must be compiled with DO_HYDRO=TRUE
#endif

#if (GRAV!=FIXED_GRAV)
#error This problem must be compiled with GRAV=FIXED_GRAV
#endif

#if (GEOM!=CARTESIAN)
#error This problem must be compiled with GEOM=CARTESIAN
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
static long int rseed;

static const double opt_L     = 10.0;
static const double opt_rho_h = 3.0;
static const double opt_rho_l = 1.0;
static const double opt_A     = 0.5;
static const double opt_g     = 2.0;
static const double opt_gamma = 5.0/3.0;
// static const double opt_norm  = 5.0e-6;
static const double opt_norm  = 3.0e-4;

static inline double opt_press()
{ 
  return 2.0*M_PI*(opt_rho_h + opt_rho_l)*opt_g*opt_L;
}

static inline double opt_rho_profile(const double rho0, const double zp)
{
  const double tmp = 1.0 - (opt_gamma - 1.0)/opt_gamma *
    (rho0*opt_g*zp)/(opt_press());
  return rho0*pow(tmp, 1.0/(opt_gamma-1.0));
}

static inline double opt_press_profile(const double rho0, const double zp)
{
  const double P0 = opt_press();
  const double rho = opt_rho_profile(rho0,zp);
  return P0*pow(rho/rho0,opt_gamma);
}

static inline double my_rand()
{
  return 2.0*ran2(&rseed) - 1.0;
}

static inline int pert_kappa_min()
{
  return n1/4;
}

static inline int pert_kappa_max()
{
  return n1/2;
}

static inline double opt_random(
  int nmodes, const double norm,
  const double * kx, const double * ky, const double * Gk,
  const double * ak, const double * bk, const double * ck, const double * dk,
  const double x, const double y)
{
  double out = 0;
  for (int i = 0; i < nmodes; ++i) {
    out += Gk[i]*ak[i]*cos(kx[i]*x)*cos(ky[i]*y);
    out += Gk[i]*bk[i]*cos(kx[i]*x)*sin(ky[i]*y);
    out += Gk[i]*ck[i]*sin(kx[i]*x)*cos(ky[i]*y);
    out += Gk[i]*dk[i]*sin(kx[i]*x)*sin(ky[i]*y);
  }
  return out/norm;
}

void init_grid()
{
  startx[0] = -opt_L/2;
  startx[1] = -opt_L/2;
  startx[2] = -opt_L;

  dx[0] = opt_L/n1;
  dx[1] = opt_L/n2;
  dx[2] = 2.0*opt_L/n3;
    
  periodic[0] = 1;
  periodic[1] = 1;
  periodic[2] = 0;
}

void init_problem()
{
  int ii, jj, kk, vv;
  double x, y, z;
  double norm;
  double *kx_vec, *ky_vec, *Gk_vec;
  double *ak_vec, *bk_vec, *ck_vec, *dk_vec;
  double x_nrm, y_nrm;

  VLOOP {
    bc[vv].lo[0] = PERIODIC;
    bc[vv].hi[0] = PERIODIC;
    bc[vv].lo[1] = PERIODIC;
    bc[vv].hi[1] = PERIODIC;
    bc[vv].lo[2] = PROB;
    bc[vv].hi[2] = PROB;
  }

  int nmodes = 0;
  for (int kx = 0; kx <= pert_kappa_max(); ++kx)
  for (int ky = 0; ky <= pert_kappa_max(); ++ky) {
    int kk = kx*kx + ky*ky;
    if (kk >= pert_kappa_min() && kk <= pert_kappa_max()) {
      ++nmodes;
    }
  }

  kx_vec = malloc(nmodes*sizeof(*kx_vec));
  ky_vec = malloc(nmodes*sizeof(*ky_vec));
  Gk_vec = malloc(nmodes*sizeof(*Gk_vec));
  ak_vec = malloc(nmodes*sizeof(*ak_vec));
  bk_vec = malloc(nmodes*sizeof(*bk_vec));
  ck_vec = malloc(nmodes*sizeof(*ck_vec));
  dk_vec = malloc(nmodes*sizeof(*dk_vec));

  memset(kx_vec, 0, nmodes*sizeof(*kx_vec));
  memset(ky_vec, 0, nmodes*sizeof(*ky_vec));
  memset(Gk_vec, 0, nmodes*sizeof(*Gk_vec));
  memset(ak_vec, 0, nmodes*sizeof(*ak_vec));
  memset(bk_vec, 0, nmodes*sizeof(*bk_vec));
  memset(ck_vec, 0, nmodes*sizeof(*ck_vec));
  memset(dk_vec, 0, nmodes*sizeof(*dk_vec));

  rseed = -myrank;

  g0[0] = 0.0;
  g0[1] = 0.0;
  g0[2] = -opt_g;

  if (mpi_io_proc()) {
    int idx = 0;
    for (int kx = 0; kx <= pert_kappa_max(); ++kx)
    for (int ky = 0; ky <= pert_kappa_max(); ++ky) {
      int kk = SQR(kx) + SQR(ky);
      if (kk >= pert_kappa_min() && kk <= pert_kappa_max()) {
        kx_vec[idx] = 2*M_PI*kx/opt_L;
        ky_vec[idx] = 2*M_PI*ky/opt_L;
        Gk_vec[idx] = sqrt(2*M_PI*opt_A*sqrt((double)kk)*opt_g/opt_L);
        ak_vec[idx] = my_rand();
        bk_vec[idx] = my_rand();
        ck_vec[idx] = my_rand();
        dk_vec[idx] = my_rand();
        norm += SQR(ak_vec[idx]) + SQR(bk_vec[idx]) +
                SQR(ck_vec[idx]) + SQR(dk_vec[idx]);
        ++idx;
      }
    }
    norm = sqrt(norm)/(opt_norm*opt_L);
  }
  MPI_Bcast(&norm,       1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(kx_vec, nmodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ky_vec, nmodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(Gk_vec, nmodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ak_vec, nmodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(bk_vec, nmodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ck_vec, nmodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(dk_vec, nmodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  gam = opt_gamma;

  ZLOOP {
    x = startx[0] + (ii + 0.5)*dx[0];
    y = startx[1] + (jj + 0.5)*dx[1];
    z = startx[2] + (kk + 0.5)*dx[2];

    double rho0 = (z < 0) ? opt_rho_l : opt_rho_h;

    double u3;
    if (fabs(z) < dx[2]) {
      u3 = opt_random(nmodes, norm,
                      kx_vec, ky_vec, Gk_vec,
                      ak_vec, bk_vec, ck_vec, dk_vec,
                      2.0*M_PI*x/opt_L, 2.0*M_PI*y/opt_L);
    } else {
      u3 = 0.0;
    }

    NDP_ELEM(sim.p,ii,jj,kk,RHO) = opt_rho_profile(rho0,z);
    NDP_ELEM(sim.p,ii,jj,kk,UU)  = opt_press_profile(rho0,z)/(opt_gamma - 1.0); 
    NDP_ELEM(sim.p,ii,jj,kk,U1)  = 0.0;
    NDP_ELEM(sim.p,ii,jj,kk,U2)  = 0.0;
    NDP_ELEM(sim.p,ii,jj,kk,U3)  = u3;
    NDP_ELEM(sim.p,ii,jj,kk,YE)  = (z > 0) ? 1.0 : 0.0;
  }

  free(dk_vec);
  free(ck_vec);
  free(bk_vec);
  free(ak_vec);
  free(Gk_vec);
  free(ky_vec);
  free(kx_vec);
}

void prob_bounds(int i, int j, int k, double *p)
{
  double z = startx[2] + (k + 0.5)*dx[2];
  double rho0 = (z < 0) ? opt_rho_l : opt_rho_h;

  p[RHO] = opt_rho_profile(rho0,z);
  p[UU]  = opt_press_profile(rho0,z)/(opt_gamma-1.0);
  p[U1]  = 0.0;
  p[U2]  = 0.0;
  p[U3]  = 0.0;
  p[YE]  = (z > 0) ? 1.0 : 0.0;

  return;
}

void analysis_preloop()
{
  return;
}

void analysis_inloop()
{
  return;
}

void analysis_postloop()
{
  return;
}
