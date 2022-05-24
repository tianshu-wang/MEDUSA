#include "decs.h"
#include "constants.h"

#define PINDEX(g,dcon,dcov) ((dcov) + SPACEDIM*((dcon) + SPACEDIM*(g)))


void p_to_u(double* p, double* u, const zone_geom* restrict g,int nhydro)
{
  // primitive to conserved
  #if (DO_HYDRO==TRUE)
  hydro_p_to_u(p, u, g->gcov[0],nhydro);
  #endif

  return;
}

void u_to_p(double* u, double* p, const zone_geom* restrict g,int nhydro,double rho_floor,double e_floor)
{
  // conserved to primitive
  #if (DO_HYDRO==TRUE)
  hydro_u_to_p(u, p, g->gcon,nhydro,rho_floor,e_floor);
  #endif

  return;
}


void calc_dvdx(int i, int j, int k, double gradv[SPACEDIM][SPACEDIM], double vavg[SPACEDIM], double xavg[NDIM])
{
  // calculates the gradient of the velocity in zone i,j,k
  // uses edge velocities from Riemann problem
  /* A.S.:  This function makes use of volume averages of linear profiles.  If some
   *   variable u can be approximated on the interval [x0,x1] by u(x) = u0 + A*(x-x0),  
   *   where u0 = u(x0), u1 = u(x1), and A = (u1-u0)/(x1-x0), then <u(x)> = A*<x> + B,
   *   where B = u0 - A*x0.  Here, <.> represents the physical volume average of the 
   *   coordinate x, and <x> = beta0/Gamma0. */
  int d,dd,ddd,jj,jp,njp,kk,kp,nkp,s,lev,dka;
  double frac,A,B,area_sum,area_frac;

  memset(gradv, 0, SQR(SPACEDIM)*sizeof(double));

  /* gradient in 0-direction */
  #if (GEOM==SPHERICAL || GEOM==CYLINDRICAL)
  njp = DJS(i)/DJS(i+1);
  jp = j*njp;
  area_sum = 0.0;
  for (jj=jp; jj<jp+njp; jj++) {
    nkp = DKS(i,j)/DKS(i+1,jj);
    kp = k*nkp;
    for (kk=kp; kk<kp+nkp; kk++) {
      area_sum += ND_ELEM_LINEAR(geom,i+1,jj,kk).area[0] + 1.0e-16;
    }
  }
  area_sum = 1.0/area_sum;
  for (jj=jp; jj<jp+njp; jj++) {
    nkp = DKS(i,j)/DKS(i+1,jj);
    kp = k*nkp;
    for (kk=kp; kk<kp+nkp; kk++) {
      area_frac = (ND_ELEM_LINEAR(geom,i+1,jj,kk).area[0]+1.0e-16)*area_sum;
      SLOOP { gradv[dd][0] += NDP_ELEM_LINEAR_F(sim_vedgedir0,0,i+1,jj,kk,dd)*area_frac; }
    }
  }
  SLOOP { 
    gradv[dd][0] -= NDP_ELEM_LINEAR_F(sim_vedgedir0,0,i,j,k,dd);
    gradv[dd][0] /= dx[0];
  }
  #else
  SLOOP { gradv[dd][0] = (NDP_ELEM_LINEAR_F(sim_vedgedir0,0,i+1,j,k,dd) - NDP_ELEM_LINEAR_F(sim_vedgedir0,0,i,j,k,dd))/dx[0]; }
  #endif
  A = gradv[0][0];
  B = NDP_ELEM_LINEAR_F(sim_vedgedir0,0,i,j,k,0) - A*(startx[0] + i*dx[0]);
  xavg[0] = beta0[i]/Gamma0[i];
  vavg[0] = A*xavg[0] + B;
  
  /* gradient in 1-direction */
  #if (NDIM>1)
  #if (GEOM==SPHERICAL || GEOM==CYLINDRICAL)
  /* outer 1-edge */
  if (DKS(i,j+1) < DKS(i,j)) {
    // Refinement boundary, 2 outer neighbors, shifted area index, original area, flux index equals area index
    SLOOP { gradv[dd][1] += 0.5*(NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j+1,2*k,dd) + NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j+1,2*k+1,dd)); }
  } else {
    // Coarsening boundary, 1 outer neighbor, shifted area index, half area, flux index equals original index
    // Regular boundary, 1 outer neighbor, original area index, original area, flux index equals original index
    SLOOP { gradv[dd][1] += NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j+1,k,dd); }
  }

  /* inner 1-edge */
  if (DKS(i,j-1) < DKS(i,j)) {
    // Coarsening boundary, 2 inner neighbors, original area index, half area, shifted flux index
    SLOOP { gradv[dd][1] -= 0.5*(NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,2*k,dd) + NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,2*k+1,dd)); }
  } else {
    // Refinement boundary, 1 inner neighbor, original area index, original area, flux index equals area index
    // Regular boundary, 1 inner neighbor, original area index, original area, flux index equals original index
    SLOOP { gradv[dd][1] -= NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,k,dd); }
  }
  SLOOP { gradv[dd][1] /= DJS(i)*dx[1]; }

  lev = 0;
  s = DJS(i);
  while (s >>= 1) lev++;
  A = gradv[1][1];
  if (DKS(i,j-1) < DKS(i,j)) {
    // Coarsening boundary, 2 inner neighbors, original area index, half area, shifted flux index
    B = 0.5*(NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,2*k,1) + NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,2*k+1,1)) - A*(startx[1] + j*DJS(i)*dx[1]);
  } else {
    // Refinement boundary, 1 inner neighbor, original area index, original area, flux index equals area index
    // Regular boundary, 1 inner neighbor, original area index, original area, flux index equals original index
    B = NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,k,1) - A*(startx[1] + j*DJS(i)*dx[1]);
  }
  xavg[1] = beta1s[lev][j]/Gamma1s[lev][j];
  vavg[1] = A*xavg[1] + B;
  #else
  SLOOP { gradv[dd][1] = (NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j+1,k,dd) - NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,k,dd))/dx[1]; }
  vavg[1] = 0.5*(NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j+1,k,1) + NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,k,1));
  #endif
  #endif
  
  /* gradient in 2-direction */
  #if (NDIM>2)
  SLOOP { gradv[dd][2] = (NDP_ELEM_LINEAR_F(sim_vedgedir2,2,i,j,k+1,dd) - NDP_ELEM_LINEAR_F(sim_vedgedir2,2,i,j,k,dd))/(DKS(i,j)*dx[2]); }
  vavg[2] = 0.5*(NDP_ELEM_LINEAR_F(sim_vedgedir2,2,i,j,k+1,2) + NDP_ELEM_LINEAR_F(sim_vedgedir2,2,i,j,k,2));
  xavg[2] = startx[2] + (k+0.5)*DKS(i,j)*dx[2]; 
  #endif
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void calc_gradv(double ** p, int i, int j, int k, double gradv[SPACEDIM][SPACEDIM])
#else
void calc_gradv(double NDP_PTR p, int i, int j, int k, double gradv[SPACEDIM][SPACEDIM])
#endif
{
  int d,dd,ddd;
  double vavg[SPACEDIM],xavg[NDIM];

  calc_dvdx(i,j,k,gradv,vavg,xavg);
  for(d=NDIM;d<SPACEDIM;d++) vavg[d] = NDP_ELEM_LINEAR(p,i,j,k,U1+d);

  for (d=0; d<SPACEDIM; d++) {
    for (dd=0; dd<SPACEDIM; dd++) {
      for (ddd=0; ddd<SPACEDIM; ddd++) {
        gradv[d][dd] += ND_ELEM_LINEAR(geom,i,j,k).conn[d][ddd][dd]*vavg[ddd];
      }
    }
  }

  return;
}

/* Note that since 
 *     (T^n_m);n = (T^n_m)_,n + \Gamma^n_nl T^l_m - \Gamma^l_mn T^n_l
 *   and since 
 *     \Gamma^n_nl = (ln g^1/2)_,l
 *   it follows that
 *     (T^n_m)_;n = g^-1/2 (g^1/2 T^n_m)_,n - \Gamma^l_nm T^n_l
 *   The first term on the right-hand side becomes the flux term.  The second term
 *   on the right-hand side becomes the geometric source term (with positive sign).
 *   Note that T[l][n] = T^n_l .
 */
GPU_PRAGMA(omp declare target)
void p_to_geom_src(double* restrict p, const double press, double* restrict src, const zone_geom* restrict gm,int nvars);
GPU_PRAGMA(omp end declare target)

void p_to_geom_src(double* restrict p, const double press, double* restrict src, const zone_geom* restrict gm,int nvars)
{
  int g,l,m,n,vv;
  double T[SPACEDIM][SPACEDIM];
  double prad,chi;
  //void rad_flux_tensor(double *p, const zone_geom* restrict gm, double* restrict Trad);

  // initialize all the source terms to zero
  // memset(src, 0, nvars*sizeof(double));
  for(int vv=0;vv<nvars;vv++) src[vv]=0;

  #if (GEOM == SPHERICAL || GEOM == CYLINDRICAL)

  #if (DO_HYDRO==TRUE)
  hydro_stress_tensor(p, press, T, gm);

  for (l=0; l<SPACEDIM; l++) {
    for (m=0; m<SPACEDIM; m++) {
      for (n=0; n<SPACEDIM; n++) {
        src[U1+m] += T[l][n]*gm->conn[l][m][n];
      }
    }
  }
  #endif

  #endif

  return;
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void get_geom_src_terms(double ** p)
#else
void get_geom_src_terms(double NDP_PTR p)
#endif
{
  static int firstc = 1;
  int ii,jj,kk,II;
  int print;

  //ZLOOP {
  //GPU_PRAGMA(omp target update to(p[:cell_count_all][:nvars],sim_eos[:cell_count_all][:NEOS],\
  //			          geom[:cell_count_all]))
  GPU_PRAGMA(omp target teams distribute parallel for)
  for(II=0;II<cell_count;II++) {
    p_to_geom_src(&NDP_ELEM_LINEAR_REORDERED(p,ii,jj,kk,0), NDP_ELEM_LINEAR_REORDERED(sim_eos,ii,jj,kk,PRESS), &NDP_ELEM_LINEAR_REORDERED(sim_src,ii,jj,kk,0), &ND_ELEM_LINEAR_REORDERED(geom,ii,jj,kk),nvars);
  }
  //GPU_PRAGMA(omp target update from(sim_src[:cell_count_all][:nvars]))

  return;
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void p_to_phys_src_start(double ** p)
#else
void p_to_phys_src_start(double NDP_PTR p)
#endif
{
  #if (GRAV!=NO_GRAV)
  gravity_start(p);
  #endif

  return;
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void p_to_phys_src_finish(double ** p, double dtpush)
#else
void p_to_phys_src_finish(double NDP_PTR p, double dtpush)
#endif
{
  int ii,jj,kk,i,j,k,g;

  TIMER_START("p_to_phys_src_finish");

  #if (GRAV!=NO_GRAV)
  // AS:  This must be placed after the sections above, since it computes a new lapse
  gravity_finish(p);
  #endif

  #if (USE_EXT_SRC==TRUE)
  ext_src(p, dtpush);
  #endif

  #if (USE_AUX_SRC==TRUE)
  aux_src(p, dtpush);
  #endif

  TIMER_STOP;

  return;
}
