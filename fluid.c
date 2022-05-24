#include "decs.h"

/* Convert primitive to conserved variables */
void hydro_p_to_u(const double* restrict p, double* restrict u, const double* restrict gcov,int nhydro)
{
  int dd;
  double vcov[SPACEDIM],vsq;

  geom_lower(&p[U1], vcov, gcov);

  u[RHO] = p[RHO];
  vsq = 0.0;
  for (dd=0; dd<SPACEDIM; dd++) {
    vsq += p[dd+U1]*vcov[dd];
    u[dd+U1] = u[RHO]*vcov[dd];
  }
  u[ETOT] = 0.5*u[RHO]*vsq + p[UU];

  for (dd=U1+SPACEDIM; dd<nhydro; dd++) {
    u[dd] = u[RHO]*p[dd];
  }

  return;
}

/* Convert conserved to primitive variables */
void hydro_u_to_p(double* restrict u, double* restrict p, const double* restrict gcon,int nhydro,double rho_floor,double e_floor)
{
  int dd;
  double vsq;

  p[RHO] = u[RHO];
  for (dd=U1; dd<U1+SPACEDIM; dd++) {
    u[dd] /= u[RHO];
  }
  geom_raise(&u[U1], &p[U1], gcon);
  vsq = 0.0;
  for (dd=0; dd<SPACEDIM; dd++) {
    vsq += u[U1+dd]*p[U1+dd];
  }
  p[UU] = u[ETOT] - 0.5*p[RHO]*vsq;
  
  #if (ENFORCE_FLOORS==TRUE)
  p[RHO] = MAX(p[RHO],rho_floor);
  p[UU ] = MAX(p[UU ],  e_floor);
  #endif

  for (dd=U1+SPACEDIM; dd<nhydro; dd++) {
    p[dd] = u[dd]/u[RHO];
  }

  return;
}


void hydro_stress_tensor(double *p, double press, double T[SPACEDIM][SPACEDIM], const zone_geom *g)
{
  // stress tensor T^i_j.  Note the location of indices
  // A.S.:  Should this be T[l][m] = T^m_l since we use vcov[l] = v_l ?
  int l,m;
  double vcov[SPACEDIM];

  geom_lower(&p[U1], vcov, g->gcov[0]);   // assumed to be at grid zone center

  for (l=0; l<SPACEDIM; l++) {
    for (m=0; m<SPACEDIM; m++) {
      T[l][m] = p[RHO]*p[U1+m]*vcov[l];
    }
    T[l][l] += press;
  }

  return;
}

void set_shock_flag(int cycle)
{
  #if (NDIM>1)
  //TIMER_START("set_shock_flag");
  // tags grid aligned shocks so HLLE can be used in directions transverse to the shock normal
  int ii,jj,kk=0,dd,grid_aligned,j,jp,njp,k,kp,nkp,dka;
  double divv,norm_dv,max_dv,dv[NDIM],shock_strength,sound_time[NDIM],vol;
  if (cycle==0) {
    //ZLOOP { 
      GPU_PRAGMA(omp target teams distribute parallel for)
      for(int II=0;II<cell_count;II++){
              GET_IJK_FROM_I(II,ii,jj,kk);
              ND_ELEM_LINEAR(sim_shock_flag,ii,jj,kk) = -1;
      }
  } 
  else {
    //ZLOOP {
      GPU_PRAGMA(omp target teams distribute parallel for firstprivate(dv,sound_time))
      for(int II=0;II<cell_count;II++){
        // calculate volume averaged velocity divergence <div(v)> in the usual finite volume fashion
        // <div(v)> = \Delta(v^i A_i)/V
        vol = ND_ELEM_LINEAR(geom,ii,jj,kk).volume;

        /* Flux differences in 0-direction */
        njp = DJS(ii)/DJS(ii+1);
        jp = jj*njp;
        dv[0] = 0.0;
        for (j=jp; j<jp+njp; j++) {
          nkp = DKS(ii,jj)/DKS(ii+1,j);
          kp = kk*nkp;
          for (k=kp; k<kp+nkp; k++) {        
            dv[0] += ND_ELEM_LINEAR(geom,ii+1,j,k).area[0]*NDP_ELEM_LINEAR_F(sim_vedgedir0,0,ii+1,j,k,0);
          }
        }
        dv[0] -= ND_ELEM_LINEAR(geom,ii,jj,kk).area[0]*NDP_ELEM_LINEAR_F(sim_vedgedir0,0,ii,jj,kk,0);
        dv[0] /= vol;
        max_dv = fabs(dv[0]);

        /* Flux differences in 1-direction */
        /* Outer 1-face */
        dv[1] = 0.0;
        if (DKS(ii,jj+1) > DKS(ii,jj)) {
          // Coarsening boundary, 1 outer neighbor, shifted area index, half area, flux index equals original index
          dv[1] += 0.5*ND_ELEM_LINEAR(geom,ii,jj+1,kk/2).area[1]*NDP_ELEM_LINEAR_F(sim_vedgedir1,1,ii,jj+1,kk,1);
        } else if (DKS(ii,jj+1) < DKS(ii,jj)) {
          // Refinement boundary, 2 outer neighbors, shifted area index, original area, flux index equals area index
          dv[1] += ND_ELEM_LINEAR(geom,ii,jj+1,2*kk  ).area[1]*NDP_ELEM_LINEAR_F(sim_vedgedir1,1,ii,jj+1,2*kk  ,1);
          dv[1] += ND_ELEM_LINEAR(geom,ii,jj+1,2*kk+1).area[1]*NDP_ELEM_LINEAR_F(sim_vedgedir1,1,ii,jj+1,2*kk+1,1);
        } else {
          // Regular boundary, 1 outer neighbor, original area index, original area, flux index equals original index
          dv[1] += ND_ELEM_LINEAR(geom,ii,jj+1,kk).area[1]*NDP_ELEM_LINEAR_F(sim_vedgedir1,1,ii,jj+1,kk,1);
        }
        
        /* Inner 1-face */
        if (DKS(ii,jj-1) < DKS(ii,jj)) {
          // Coarsening boundary, 2 inner neighbors, original area index, half area, shifted flux index
          dv[1] -= 0.5*ND_ELEM_LINEAR(geom,ii,jj,kk).area[1]*NDP_ELEM_LINEAR_F(sim_vedgedir1,1,ii,jj,2*kk  ,1);
          dv[1] -= 0.5*ND_ELEM_LINEAR(geom,ii,jj,kk).area[1]*NDP_ELEM_LINEAR_F(sim_vedgedir1,1,ii,jj,2*kk+1,1);
        } else {
          // Refinement boundary, 1 inner neighbor, original area index, original area, flux index equals area index
          // Regular boundary, 1 inner neighbor, original area index, original area, flux index equals original index
          dv[1] -= ND_ELEM_LINEAR(geom,ii,jj,kk).area[1]*NDP_ELEM_LINEAR_F(sim_vedgedir1,1,ii,jj,kk,1);
        }
        dv[1] /= vol;
        max_dv = MAX(max_dv,fabs(dv[1]));

        #if (NDIM>2)
        /* Flux differences in 2-direction */
        dv[2] = (ND_ELEM_LINEAR(geom,ii,jj,kk+1).area[2]*NDP_ELEM_LINEAR_F(sim_vedgedir2,2,ii,jj,kk+1,2)
              -  ND_ELEM_LINEAR(geom,ii,jj,kk  ).area[2]*NDP_ELEM_LINEAR_F(sim_vedgedir2,2,ii,jj,kk  ,2))/vol;
        max_dv = MAX(max_dv,fabs(dv[2]));
        #endif

        // normalize to unit magnitude and calculate sound crossing time
        divv = 0.0;
        norm_dv = 0.0;
        DLOOP {
          divv += dv[dd];
          norm_dv += dv[dd]*dv[dd];
          sound_time[dd] = dx[dd]*ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][dd]/NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,CS);
          #if (NDIM>1 && (GEOM==SPHERICAL || GEOM==CYLINDRICAL))
          if (dd==1) sound_time[dd] *= DJS(ii);
          #endif
          #if (NDIM>2 && GEOM==SPHERICAL)
          if (dd==2) sound_time[dd] *= DKS(ii,jj);
          #endif
        }
        norm_dv = sqrt(norm_dv);

        // check if velocity is converging and the convergence is dominated by a single direction
        grid_aligned = -1;
        DLOOP {
          dv[dd] /= -(norm_dv + 1.0e-10);
          if (dv[dd] > 0.9) {
            grid_aligned = dd;
            break;
          }
        }

        // check if the convergence rate is trans/supersonic.  be liberal here.
        ND_ELEM_LINEAR(sim_shock_flag,ii,jj,kk) = -1;
        if (grid_aligned >= 0) { 
          shock_strength = -divv * sound_time[grid_aligned];
          if (shock_strength > 0.5) {
            ND_ELEM_LINEAR(sim_shock_flag,ii,jj,kk) = grid_aligned;
          }
        }
      }
  }
  //TIMER_STOP;
  #endif /* NDIM>1 */
  return;
}

int transverse_shock(int ND_PTR ijk_to_I,int *sim_shock_flag, int *dj, int **dk,int i, int j, int k, int dir)
{
  int km,nkm,kk,orflag;
  
  // check shock_flag to see if there is a grid aligned shock orthogonal to dir
  // WARNING: we're only dealing with shocks along direction 0.  We need more logic for the other directions
  // A.S.:  shock_flag=={-1,0} if {no,yes} shock
  if (dir == 1) {
    nkm = DKS(i,j)/DKS(i,j-1);
    km  = k*nkm;
    orflag = (ND_ELEM_LINEAR(sim_shock_flag,i,j,k) == 0);
    for (kk=km; kk<km+nkm; kk++) {
      orflag = (orflag || ND_ELEM_LINEAR(sim_shock_flag,i,j-1,kk) == 0);
    }
    return orflag;
    // if (ND_ELEM_LINEAR(sim_shock_flag,i,j,k) == 0 || ND_ELEM_LINEAR(sim_shock_flag,i,j-1,km) == 0) return 1;
  } else if (dir == 2) {
    if (ND_ELEM_LINEAR(sim_shock_flag,i,j,k) == 0 || ND_ELEM_LINEAR(sim_shock_flag,i,j,k-1) == 0) return 1;
  }

  return 0;
}
