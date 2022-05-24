#include "../decs.h"

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

void init_grid()
{
  startx[0] = 0.0;
  startx[1] = 0.0;
  startx[2] = 0.0;
  dx[0] = 0.5/n1;
  dx[1] = 1.0/n2;
  dx[2] = 0.5/n3;

  periodic[0] = 1;
  periodic[1] = 0;
  periodic[2] = 1;

  return;
}

void init_problem()
{
  double x,y,ybound;
  double rhotop,rhobot;
  double pert_amp,hsmooth;
  double pbase;
  double xl,xh,yl,yh;
  int ii,jj,kk,nsub,vv,i,j;

  gam      = 1.4;
  rhotop   = 2.0;
  rhobot   = 1.0;
  pert_amp = 0.01;
  hsmooth  = 0.005;
  pbase    = 10.0/7.0 + 1.0/4.0;
  nsub     = 50;

  g0[0] = 0.0;
  g0[1] = -0.5;
  g0[2] = 0.0;

  VLOOP {
    bc[vv].lo[0] = PERIODIC;
    bc[vv].hi[0] = PERIODIC;
    bc[vv].lo[1] = REFLECT;
    bc[vv].hi[1] = REFLECT;
    #if NDIM==3
    bc[vv].lo[2] = PERIODIC;
    bc[vv].hi[2] = PERIODIC;
    #endif
  }

  ZLOOP {
    x = startx[0] + (ii+0.5)*dx[0];
    y = startx[1] + (jj+0.5)*dx[1];

    xl = startx[0] + ii*dx[0];
    xh = xl + dx[0];
    yl = startx[0] + jj*dx[1];
    yh = yl + dx[1];

    if (yl > 0.5+pert_amp || yh < 0.5-pert_amp) {
      if (y < 0.5) {
        NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhobot;
        NDP_ELEM(sim.p,ii,jj,kk,UU ) = (pbase + rhobot*g0[1]*y)/(gam-1.0);
      } else {
        NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhotop;
        NDP_ELEM(sim.p,ii,jj,kk,UU ) = (pbase + rhobot*g0[1]*0.5 + rhotop*g0[1]*(y-0.5))/(gam-1.0);
      }
    } else {
       
      NDP_ELEM(sim.p,ii,jj,kk,RHO) = 0.0;
      NDP_ELEM(sim.p,ii,jj,kk,UU ) = 0.0;
      for (i=0; i<nsub; i++) {
        x = xl + (i+0.5)*dx[0]/nsub;
        ybound = 0.5*pert_amp*(cos(8.0*M_PI*x)) + 0.5;
        for (j=0;j<nsub;j++) {
          y = yl + (j+0.5)*dx[1]/nsub;
          if (y < ybound) {
            NDP_ELEM(sim.p,ii,jj,kk,RHO) += rhobot + 0.5*(rhotop-rhobot) * (1.0 + tanh((y-ybound)/hsmooth));
            NDP_ELEM(sim.p,ii,jj,kk,UU ) += (pbase + rhobot*g0[1]*y)/(gam-1.0);
          } else {
            NDP_ELEM(sim.p,ii,jj,kk,RHO) += rhobot + 0.5*(rhotop-rhobot) * (1.0 + tanh((y-ybound)/hsmooth));
            NDP_ELEM(sim.p,ii,jj,kk,UU ) += (pbase + rhobot*g0[1]*ybound + rhotop*g0[1]*(y-ybound))/(gam-1.0);
          }
        }
	
      }
      NDP_ELEM(sim.p,ii,jj,kk,RHO) /= SQR(nsub);
      NDP_ELEM(sim.p,ii,jj,kk,UU ) /= SQR(nsub);
    }

    NDP_ELEM(sim.p,ii,jj,kk,YE) = NDP_ELEM(sim.p,ii,jj,kk,RHO) - 1.0;
    NDP_ELEM(sim.p,ii,jj,kk,U1) = 0.0;
    NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.0;
    NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.0;

    /*        
    if (y > 0.5) NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhotop;
    else NDP_ELEM(sim.p,ii,jj,kk,RHO) = rhobot;
    //NDP_ELEM(sim.p,ii,jj,kk,UU) = (10./7. - g0[1]*NDP_ELEM(sim.p,ii,jj,kk,RHO)*(y - 0.5))/(gam-1.);
    if (y < 0.5) NDP_ELEM(sim.p,ii,jj,kk,UU) = (2.5 + g0[1]*rhobot*y)/(gam-1);
    else NDP_ELEM(sim.p,ii,jj,kk,UU) = (2.5 + g0[1]*0.5*rhobot + g0[1]*(y-0.5)*rhotop)/(gam-1);
    NDP_ELEM(sim.p,ii,jj,kk,U1) = 0.;
    if (y > 0.3 && y < 0.7) {
    NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.1*cos(8*M_PI*(x+0.25))*(1 + cos(5*M_PI*(y-0.5)))/4./NDP_ELEM(sim.p,ii,jj,kk,RHO);
    } else {
    NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.;
    }
    NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;
    */
  }


	return;
}

void prob_bounds(int i, int j, int k, double *p)
{
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
