#include "decs.h"
#include "constants.h"

void init_tracers_by_mass(long int *rseed) {
  int ii,jj,kk,i,j;
  double cxlo[3],cxhi[3];
  double *ms = malloc_rank1(n1, sizeof(double));
  double *m_part_r = malloc_rank1(n1, sizeof(double));
  double *mth = malloc_rank1(n2, sizeof(double));
  double *m_part_th = malloc_rank1(n1, sizeof(double));
  for(ii=0;ii<n1;ii++) {
    ms[ii] = 0.;
    m_part_r[ii] = 0.;
  }
  for(jj=0;jj<n2;jj++) {
    mth[jj] = 0.;
    m_part_th[jj] = 0.;
  }
  ZLOOP {
    double dm = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO)*ND_ELEM_LINEAR(geom,ii,jj,kk).volume;
    ms[ii] += dm;
    double cthhi = cos(th_of_x(thx_info,startx[1] + jj*DJS(ii)*dx[1]));
    double cthlo;
    double denom = 1. / (cos(th_of_x(thx_info,startx[1] + (jj+1)*DJS(ii)*dx[1])) - cthhi);
    for (j=jj*DJS(ii); j<(jj+1)*DJS(ii); ++j) {
      cthlo = cthhi;
      cthhi = cos(th_of_x(thx_info,startx[1] + (j+1)*dx[1]));
      mth[j] += dm * (cthhi - cthlo) * denom;
    }
  }
  mpi_global_reduce(ms, n1);
  mpi_global_reduce(mth, n2);

  int i_inner_tracers = -1;
  int i_outer_tracers = -1;
  double mass_tot = 0.;
  double mass_inside = 0.;
  double mass_outside = 0.;

  mass_outside_tracers *= MSUN;
  mass_inside_tracers *= MSUN;
  for(ii=0;ii<n1;ii++) {
      mass_tot += ms[ii];
      if(i_inner_tracers < 0 && mass_tot > mass_inside_tracers) {
          i_inner_tracers = ii;
          mass_inside = mass_tot - ms[ii];
      }
      if(i_outer_tracers < 0 && mass_tot > mass_outside_tracers){
          i_outer_tracers = ii;
          mass_outside = mass_tot;
      }
  }
  if (i_outer_tracers < 0) {
    mass_outside = mass_tot;
    i_outer_tracers = n1 - 1;
  }

  double total_traced_mass = mass_outside - mass_inside;
  double tracer_mass_target = total_traced_mass/n_tracer_target;

  n_tracer_current = 0;
  tracers = NULL;
  ZLOOP {
    if(ii>=i_inner_tracers && ii <= i_outer_tracers) {
      cxlo[0] = r_of_x(rx_info,startx[0] +  ii   *dx[0]);
      cxhi[0] = r_of_x(rx_info,startx[0] + (ii+1)*dx[0]);
      #if(NDIM>1)
      cxlo[1] = startx[1] +  jj   *DJS(ii)*dx[1];
      cxhi[1] = startx[1] + (jj+1)*DJS(ii)*dx[1];
      #if(NDIM==3)
      cxlo[2] = startx[2] +  kk   *DKS(ii,jj)*dx[2];
      cxhi[2] = startx[2] + (kk+1)*DKS(ii,jj)*dx[2];
      #endif
      #endif
      double mcell = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO)*ND_ELEM_LINEAR(geom,ii,jj,kk).volume;
      int i_in_cell = (int)(mcell/tracer_mass_target);
      double msamp = mcell - i_in_cell*tracer_mass_target;
      int n_in_cell = i_in_cell + (ran2(rseed) < msamp/tracer_mass_target ? 1 : 0);
      for(i=0;i<n_in_cell;i++) {
        tracer_t *newt = malloc_rank1(1, sizeof(tracer_t));
        newt->mass = mcell/n_in_cell;
        for(j=0;j<NDIM;j++) {
          newt->x[j] = cxlo[j] + ran2(rseed)*(cxhi[j]-cxlo[j]);
          if (j==0) {
            m_part_r[ii] += newt->mass;
          } else if (j==1) {
            int j0 = jj+(int)((newt->x[j]-cxlo[j])/dx[j]);
            //printf("th=%.3g, dth=%.3g, j=%d\n", newt->x[j], (newt->x[j]-cxlo[j])/dx[j], j0);
            m_part_th[jj+(int)((newt->x[j]-cxlo[j])/dx[j])] += newt->mass;
          }
        }
        #if ((RCOORD==UNIFORM) || (RCOORD==SINH))
        newt->x[0] = x_of_r(newt->x[0]);
        #else
        newt->x[0] = comp_x_from_r(newt->x[0], startx[0] + ii*dx[0],
                                   startx[0] + (ii+1)*dx[0]);
        #endif
        newt->id = n_tracer_current+i;
        newt->active = 1;
        newt->next = tracers;
        tracers = newt;
      }
      n_tracer_current += n_in_cell;
    }
  }
  // correct tracer masses
  mpi_global_reduce(m_part_r, n1);
  mpi_global_reduce(m_part_th, n2);
  tracer_t *ptr = tracers;
  while (ptr) {
    i = (int)((ptr->x[0] - startx[0]) / dx[0]);
    #if (NDIM > 1)
    j = (int)((ptr->x[1] - startx[1]) / dx[1]);
    ptr->mass *= (ms[i] / m_part_r[i]) * (mth[j] / m_part_th[j]);
    #else
    ptr->mass *= (ms[i] / m_part_r[i]);
    #endif
    ptr = ptr->next;
  }

  n_tracer_global = n_tracer_current;
  mpi_global_reduce_uint(&n_tracer_global, 1);
  free(ms);
  free(mth);
  free(m_part_r);
  free(m_part_th);
}
