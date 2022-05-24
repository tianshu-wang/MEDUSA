#include "../decs.h"

#if (GRAV!=NO_GRAV)
#error This problem must be compiled with GRAV=NO_GRAV
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
  dx[0] = 4.0/n1;
  dx[1] = 2.0/n2;
  dx[2] = 1.0/n3;

  periodic[0] = periodic[1] = periodic[2] = 1;

  return;
}

void init_problem()
{
  int i,j,ii,jj,kk,vv,dd;
  double x,y,t1,t2;
  double pert_amp,rho,v,w;

  gam = 5.0/3.0;

  double drho = 1.0;
  pert_amp = 0.01;
  double rho0 = 1.0;
  double u0 = 1.0;

  VLOOP {
    DLOOP {
      bc[vv].lo[dd] = PERIODIC;
      bc[vv].hi[dd] = PERIODIC;
    }
  }

  double ybound = 3.0 - (1.0 + 1.0/sqrt(2.0));
  int sub = 16;
  double y1 = 0.5;
  double y2 = 1.5;
  double a = 0.05;
  double sig = 0.2;

  ZLOOP {
    NDP_ELEM(sim.p,ii,jj,kk,RHO) = 0.0;
    NDP_ELEM(sim.p,ii,jj,kk,U1 ) = 0.0;
    NDP_ELEM(sim.p,ii,jj,kk,U2 ) = 0.0;
    for (i=0; i<sub; i++) {
      for (j=0; j<sub; j++) {
        x = startx[0] + ii*dx[0] + (i+0.5)*dx[0]/sub;
        y = startx[1] + jj*dx[1] + (j+0.5)*dx[1]/sub;
        t1 = tanh((y-y1)/a);
        t2 = tanh((y-y2)/a);
        rho = 1.0 + drho/rho0 * 0.5 * (t1 - t2);
        v = u0*(t1 - t2 - 1.0);
        w = pert_amp*sin(M_PI*x)*(exp(-SQR(y-y1)/SQR(sig)) + exp(-SQR(y-y2)/SQR(sig)));
        NDP_ELEM(sim.p,ii,jj,kk,RHO) += rho;
        NDP_ELEM(sim.p,ii,jj,kk,U1 ) += v;
        NDP_ELEM(sim.p,ii,jj,kk,U2 ) += w;
      }
    }
    NDP_ELEM(sim.p,ii,jj,kk,RHO) /= SQR(sub);
    NDP_ELEM(sim.p,ii,jj,kk,U1 ) /= SQR(sub);
    NDP_ELEM(sim.p,ii,jj,kk,U2 ) /= SQR(sub);
    NDP_ELEM(sim.p,ii,jj,kk,UU )  = 10.0/(gam-1.0);
  }

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
