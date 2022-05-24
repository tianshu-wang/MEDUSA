#include "../decs.h"

#if (DO_HYDRO!=TRUE)
#error This problem must be compiled with DO_HYDRO=TRUE
#endif

#if (GRAV!=NO_GRAV)
#error This problem must be compiled with GRAV=NO_GRAV
#endif

#if (DO_RADIATION!=FALSE)
#error This problem must be compiled with DO_RADIATION=FALSE
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

void init_grid()
{
  startx[0] = 0.0;
  dx[0] = 1.0/n1;

  return;
}

void init_problem()
{
  int ii,jj,kk,dd,vv;
  double xl,xr;
  double rhol,el,ul;
  double rhor,er,ur;

  gam = 1.4;

  rhol = 1.0;
  el = 1.0/(gam-1.0);
  ul = 0.0;

  rhor = 0.125;
  er = 0.1/(gam-1.0);
  ur = 0.0;

  fprintf(stderr,"[init_problem]:  my_grid_dims = %d\n", my_grid_dims[0]);

  // set the boundary conditions
  VLOOP {
    bc[vv].lo[0] = OUTFLOW;
    bc[vv].hi[0] = OUTFLOW;
  }


  ZLOOP {
    xl = ii*dx[0];
    xr = xl+dx[0];

    if (xr < 0.5) {
      NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhol;
      NDP_ELEM(sim.p,ii,jj,kk,UU ) = el;
      DLOOP {
        NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = ul;
      }
    } else if (xl < 0.5) {
      NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhol*(0.5-xl)/dx[0] + rhor*(xr-0.5)/dx[0];
      NDP_ELEM(sim.p,ii,jj,kk,UU ) = el*(0.5-xl)/dx[0] + er*(xr-0.5)/dx[0];
      DLOOP {
        NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = ul*(0.5-xl)/dx[0] + ur*(xr-0.5)/dx[0];
      }
    } else {
      NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhor;
      NDP_ELEM(sim.p,ii,jj,kk,UU ) = er;
      DLOOP {
        NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = ur;
      }
    }

    for (dd=NDIM; dd<SPACEDIM; dd++) {
      NDP_ELEM(sim.p,ii,jj,kk,U1+dd) = 0.0;
    }

  }

  reset_boundaries(sim.p);
  update_eos(sim.p, sim.eos);

  return;
}

void prob_bounds(int i, int j, int k, double *p)
{
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
