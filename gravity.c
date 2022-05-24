#include "decs.h"
#include "constants.h"

// #define THIRD_ORDER_MONOPOLE

#define LMAX1 ((multipole_lmax+1))
#define LM_INDEX(l,m)       ((l)*((l)+1) + (m))
#define ILM_INDEX(i,l,m)    ((i)*LMAX1*LMAX1 + LM_INDEX((l),(m)))

#if (GRAV==USER_GRAV)
/* This function calls a user-defined grav_accel() function, to be defined in
 * in the problem initialization file, that must return the gravitational
 * acceleration as a *contravariant* vector.  For the total energy source term,
 * the components must first be lowered before being dotted with the contravariant
 * velocity in the primitive variables.  For the momentum source term, the conserved
 * momenta are already covariant, so must add gravity source term as covariant as well.
 */
#if (USE_LINEAR_ALLOCATION==TRUE)
void user_gravity(double ** p)
#else
void user_gravity(double NDP_PTR p)
#endif
{
  int ii,jj,kk,dd;

  ZLOOP {
    grav_accel(ii,jj,kk,g0);
    NDP_ELEM_LINEAR(sim_src,ii,jj,kk,ETOT) += NDP_ELEM_LINEAR(p,ii,jj,kk,RHO)*geom_dot(&NDP_ELEM_LINEAR(p,ii,jj,kk,U1),g0,ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0]);
    DLOOP {
      NDP_ELEM_LINEAR(sim_src,ii,jj,kk,P1+dd) += NDP_ELEM_LINEAR(p,ii,jj,kk,RHO)*g0[dd]*ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][dd];
    }
  }

  return;
}
#endif

#if (GRAV==FIXED_GRAV)
/* This function uses a time-constant gravitational acceleration defined in the problem
 * initialization file as a *contravariant* vector.  For the total energy source term,
 * the components must first be lowered before being dotted with the contravariant
 * velocity in the primitive variables.  For the momentum source term, the conserved
 * momenta are already covariant, so must add gravity source term as covariant as well.
 */
#if (USE_LINEAR_ALLOCATION==TRUE)
void fixed_gravity(double ** p)
#else
void fixed_gravity(double NDP_PTR p)
#endif
{
  int ii,jj,kk,dd;

  ZLOOP {
    NDP_ELEM_LINEAR(sim_src,ii,jj,kk,ETOT) += NDP_ELEM_LINEAR(p,ii,jj,kk,RHO)*geom_dot(&NDP_ELEM_LINEAR(p,ii,jj,kk,U1),g0,ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0]);
    DLOOP {
      NDP_ELEM_LINEAR(sim_src,ii,jj,kk,P1+dd) += NDP_ELEM_LINEAR(p,ii,jj,kk,RHO)*g0[dd]*ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][dd];
    }
  }

  return;
}
#endif

#if (GRAV==PRESCRIBED_GRAV)
#if (USE_LINEAR_ALLOCATION==TRUE)
void prescribed_gravity(double ** p)
#else
void prescribed_gravity(double NDP_PTR p)
#endif
{
  int ii,jj,kk,dd;
  double xm,r,rhat[SPACEDIM],gr,g[SPACEDIM],v[SPACEDIM];
  double vl,vr=0;
  //double ijk_to_r(int i, int j, int k, double rhat[]);

  ZLOOP {
    /*r = ijk_to_r(ii,jj,kk,rhat);
    gr = -GNEWT*M_prescribed/(r*r);
    SLOOP {
    g[dd] = gr*rhat[dd];
    v[dd] = NDP_ELEM_LINEAR(p,ii,jj,kk,U1+dd);
    }
    NDP_ELEM_LINEAR(sim_src,ii,jj,kk,ETOT) += NDP_ELEM_LINEAR(p,ii,jj,kk,RHO)*DOT(v, g);
    DLOOP {
    NDP_ELEM_LINEAR(sim_src,ii,jj,kk,U1+dd) += NDP_ELEM_LINEAR(p,ii,jj,kk,RHO)*g[dd];
    }*/

    //vl = NDP_ELEM(sim_vedge[0],ii,jj,kk,0);
    vl = NDP_ELEM_LINEAR_F(sim_fdir0,0,ii,jj,kk,0)/ND_ELEM_LINEAR(geom,ii,jj,kk).area[0];  // AS:  likely to cause a problem at origin...
    int jp = jj;
    int kp = kk;
    if (DJS(ii+1) != DJS(ii)) {
      jp = istart[1] + jj*2;  // BROKEN
      #if (GEOM==SPHERICAL && NDIM==3)
      kp = istart[2] + kk*2;  // BROKEN
      //vr = 0.25* (NDP_ELEM(sim_vedge[0],ii+1,jp,kp,0) + NDP_ELEM(sim_vedge[0],ii+1,jp+1,kp,0)
      //            + NDP_ELEM(sim_vedge[0],ii+1,jp,kp+1,0) + NDP_ELEM(sim_vedge[0],ii+1,jp+1,kp+1,0));
      vr = 0.25* (NDP_ELEM_LINEAR_F(sim_fdir0,0,ii+1,jp,kp,0)/ND_ELEM_LINEAR(geom,ii+1,jp,kp).area[0] + NDP_ELEM_LINEAR_F(sim_fdir0,0,ii+1,jp+1,kp,0)/ND_ELEM_LINEAR(geom,ii+1,jp+1,kp).area[0]
        + NDP_ELEM_LINEAR_F(sim_fdir0,0,ii+1,jp,kp+1,0)/ND_ELEM_LINEAR(geom,ii+1,jp,kp+1).area[0] + NDP_ELEM_LINEAR_F(sim_fdir0,0,ii+1,jp+1,kp+1,0)/ND_ELEM_LINEAR(geom,ii+1,jp+1,kp+1).area[0]);
      #elif (NDIM==2)
      //vr = 0.5*(NDP_ELEM(sim_vedge[0],ii+1,jp,kp,0) + NDP_ELEM(sim_vedge[0],ii+1,jp+1,kp,0));
      vr = 0.5*(NDP_ELEM_LINEAR_F(sim_fdir0,0,ii+1,jp,kp,0)/ND_ELEM_LINEAR(geom,ii+1,jp,kp).area[0] + NDP_ELEM_LINEAR_F(sim_fdir0,0,ii+1,jp+1,kp,0)/ND_ELEM_LINEAR(geom,ii,jp+1,kp).area[0]);
      #endif
    } else {
      //vr = NDP_ELEM(sim_vedge[0],ii+1,jp,kp,0);
      vr = NDP_ELEM_LINEAR_F(sim_fdir0,0,ii+1,jp,kp,0)/ND_ELEM_LINEAR(geom,ii+1,jp,kp).area[0];
    }

    double xl = startx[0] + ii*dx[0];
    double x0 = beta0[ii]/Gamma0[ii];

    //double gl = -GNEWT*M_prescribed/pow(r_of_x(rx_info,xl),2) * dr_dx(xl);
    //double gr = -GNEWT*M_prescribed/pow(r_of_x(rx_info,xr),2) * dr_dx(xr);
    #if (PRESCRIBED_GRAV==TRUE && PN_POTENTIAL==TRUE)
    double g0 = -GNEWT*M_prescribed/pow(r_of_x(rx_info,x0) - Rschw,2) * dr_dx(x0);
    #else
    double g0 = -GNEWT*M_prescribed/pow(r_of_x(rx_info,x0),2) * dr_dx(x0);
    #endif

    //double A = (gr - gl)/dx[0];
    //double B = gl - A*xl;
    double momsrc = NDP_ELEM_LINEAR(p,ii,jj,kk,RHO)*g0;//(A*beta0[ii]/Gamma0[ii] + B);
    NDP_ELEM_LINEAR(sim_src,ii,jj,kk,U1) += momsrc;

    double A = (vr - vl)/dx[0];
    double B = vl - A*xl;
    double v1 = A*beta0[ii]/Gamma0[ii] + B;

    //NDP_ELEM_LINEAR(sim_src,ii,jj,kk,ETOT) += v1*momsrc;
    NDP_ELEM_LINEAR(sim_src,ii,jj,kk,ETOT) += v1*g0;

  }
  //exit(1);

  return;
}
#endif

#if (GEOM==SPHERICAL)
static double ms[NSIZE_GRID],*mr,*grav,*dgrav;
#if (GR_MONOPOLE==TRUE)
static double rhoavg[NSIZE_GRID],pgasavg[NSIZE_GRID],pradavg[NSIZE_GRID],uavg[NSIZE_GRID],vavg[NSIZE_GRID],eradavg[NSIZE_GRID],vdF[NSIZE_GRID],*gtov;
#endif
#if (USE_LINEAR_ALLOCATION==TRUE)
void monopole_gravity_start(double ** p)
#else
void monopole_gravity_start(double NDP_PTR p)
#endif
{
  // This routine computes the gravitational momentum and energy source terms with a monopole approximation
  // Assumes that the first dimension is radial.

  static int firstc = 1;

  //TIMER_START("monopole_gravity_start");

  int ii,jj,kk,dd,g;
  double momsrc,Etot,vdotF,Fmag,fred,chi,Prr,vol;
  double Ftot[SPACEDIM];

  if (firstc) {
    //ms            = malloc_rank1(n1,   sizeof *ms           );
    mr            = malloc_rank1(n1+1, sizeof *mr           );
    grav          = malloc_rank1(n1+1, sizeof *grav         );
    dgrav         = malloc_rank1(n1,   sizeof *dgrav        );
    #if (GR_MONOPOLE==TRUE)
    //rhoavg        = malloc_rank1(n1,   sizeof *rhoavg       );
    //pgasavg       = malloc_rank1(n1,   sizeof *pgasavg      );
    //pradavg       = malloc_rank1(n1,   sizeof *pradavg      );
    //uavg          = malloc_rank1(n1,   sizeof *uavg         );
    //vavg          = malloc_rank1(n1,   sizeof *vavg         );
    //eradavg       = malloc_rank1(n1,   sizeof *eradavg      );
    //vdF           = malloc_rank1(n1,   sizeof *vdF          );
    gtov          = malloc_rank1(n1+1, sizeof *gtov         );
    gr_grav       = malloc_rank1(n1,   sizeof *gr_grav      );
    #endif
    firstc=0;
  }

  // assumes that the first dimension is radial
  for (ii=0; ii<NSIZE_GRID; ii++) {
    ms[ii]      = 0.0;
    #if (GR_MONOPOLE==TRUE)
    rhoavg[ii]  = 0.0;
    pgasavg[ii] = 0.0;
    pradavg[ii] = 0.0;
    uavg[ii]    = 0.0;
    vavg[ii]    = 0.0;
    eradavg[ii] = 0.0;
    vdF[ii]     = 0.0;
    #endif
  }

  #if (GR_MONOPOLE==TRUE)
  //ZLOOP {
  GPU_PRAGMA(omp target data map(tofrom:rhoavg,pgasavg,pradavg,uavg,\
			                vavg,eradavg,vdF)){
  GPU_PRAGMA(omp target teams distribute parallel for
	     //reduction(+:rhoavg[:NSIZE_GRID],pgasavg[:NSIZE_GRID],pradavg[:NSIZE_GRID],\
	     //          uavg[:NSIZE_GRID],vavg[:NSIZE_GRID],eradavg[:NSIZE_GRID],vdF[:NSIZE_GRID])
	     )
    ZLOOP { //only parallel over ii, so no need to reduce the arrays. 
    //for(int II=0;II<cell_count;II++) {
      //GET_IJK_FROM_I(II,ii,jj,kk);
      #if (DO_RADIATION==TRUE)
      Etot = 0.0;
      SLOOP { Ftot[dd] = 0.0; }
      GLOOP {
        Etot += NDP_ELEM_LINEAR(p,ii,jj,kk,irad1+g);
        SLOOP { Ftot[dd] += NDP_ELEM_LINEAR(p,ii,jj,kk,ifrad1+g*NDIM+dd); }
      }
      Etot /= SQR(CLIGHT);
      vdotF = geom_dot(Ftot, &NDP_ELEM_LINEAR(p,ii,jj,kk,U1), ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0])/pow(CLIGHT,4);
      // calculate Prr for the M1 closure
      Fmag  = sqrt(geom_dot(Ftot,Ftot,ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0]));
      fred  = Fmag/(CLIGHT*fabs(Etot) + 1.0e-16);
      fred  = MAX(0.0,MIN(1.0,fred));
      chi   = (3.0 + 4.0*SQR(fred))/(5.0 + 2.0*sqrt(4.0 - 3.0*SQR(fred)));
      Prr   = 0.5*Etot*((1.0-chi) + (3.0*chi-1.0)*SQR(Ftot[0])*ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][0]/(SQR(Fmag) + 1.0e-16));
      #else
      Etot  = 0.0;
      vdotF = 0.0;
      Prr   = 0.0;
      #endif
      vol          = ND_ELEM_LINEAR(geom,ii,jj,kk).volume;
      rhoavg[ii]  += vol*NDP_ELEM_LINEAR(p,ii,jj,kk,RHO);
      pgasavg[ii] += vol*NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,PRESS);
      pradavg[ii] += vol*Prr;
      uavg[ii]    += vol*NDP_ELEM_LINEAR(p,ii,jj,kk,UU);
      vavg[ii]    += vol*NDP_ELEM_LINEAR(p,ii,jj,kk,U1);
      eradavg[ii] += vol*Etot;
      vdF[ii]     += vol*vdotF;
    }
  }
    #else  /* NOT GR_MONOPOLE */
    //ZLOOP { 
    GPU_PRAGMA(omp target data map(tofrom:ms)){
    GPU_PRAGMA(omp target teams distribute parallel for reduction(+:ms[:NSIZE_GRID]))
    for(int II=0;II<cell_count;II++) {
      GET_IJK_FROM_I(II,ii,jj,kk);
      ms[ii] += ND_ELEM_LINEAR(geom,ii,jj,kk).volume*NDP_ELEM_LINEAR(p,ii,jj,kk,RHO); 
    }
  }
  #endif  /* ifdef GR_MONOPOLE */

  #if (GR_MONOPOLE==TRUE)
  mpi_global_ireduce_start(rhoavg,  n1, 0);
  mpi_global_ireduce_start(pgasavg, n1, 1);
  mpi_global_ireduce_start(pradavg, n1, 2);
  mpi_global_ireduce_start(uavg,    n1, 3);
  mpi_global_ireduce_start(vavg,    n1, 4);
  mpi_global_ireduce_start(eradavg, n1, 5);
  mpi_global_ireduce_start(vdF,     n1, 6);
  #else
  mpi_global_ireduce_start(ms,      n1, 0);
  #endif

  //TIMER_STOP;
  return;
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void monopole_gravity_finish(double ** p)
#else
void monopole_gravity_finish(double NDP_PTR p)
#endif
{
  static int firstc = 1;
  int ii,jj=0,kk=0,g,iter,j,k,jp,kp,njp,nkp;
  double momsrc,gtov_avg,mlast,mstar;
  double Vshell,xl,A,B,vedge,vr,vl,area;
  double pgasedge,pradedge,rhoedge,uedge;
  double C,xr;

  TIMER_START("monpole_gravity_finish");

  #if (GR_MONOPOLE==TRUE)
  mpi_global_ireduce_finish(rhoavg,  n1, 0);
  mpi_global_ireduce_finish(pgasavg, n1, 1);
  mpi_global_ireduce_finish(pradavg, n1, 2);
  mpi_global_ireduce_finish(uavg,    n1, 3);
  mpi_global_ireduce_finish(vavg,    n1, 4);
  mpi_global_ireduce_finish(eradavg, n1, 5);
  mpi_global_ireduce_finish(vdF,     n1, 6);

  // need to get volume averages
  for (ii=0; ii<n1; ii++) {
    Vshell       = 1.0/(4.0*M_PI/3.0 * (pow(redge[ii+1],3) - pow(redge[ii],3)));
    rhoavg[ii]  *= Vshell;
    pgasavg[ii] *= Vshell;
    pradavg[ii] *= Vshell;
    uavg[ii]    *= Vshell;
    vavg[ii]    *= Vshell;
    eradavg[ii] *= Vshell;
    vdF[ii]     *= Vshell;
  }

  // now iterate to converge on Mtov and Gamma
  mr[0]   = 0.0;
  grav[0] = 0.0;
  if (firstc) {
    // get a reasonable estimate for gtov
    gtov[0] = 1.0;
    mstar   = 0.0;
    for (ii=1; ii<=n1; ii++) {
      Vshell   = 4.0*M_PI/3.0 * (pow(redge[ii],3) - pow(redge[ii-1],3));
      mstar   += Vshell*rhoavg[ii-1];
      gtov[ii] = sqrt(1.0 - 2.0*GNEWT*mstar/(SQR(CLIGHT)*redge[ii]));
    }
    firstc = 0;
  }

  iter = 0;
  while (1) {
    mstar = 0.0;
    for (ii=1; ii<=n1; ii++) {
      xl = startx[0] + ii*dx[0];
      A  = (gtov[ii] - gtov[ii-1])/dx[0];
      B  = gtov[ii-1] - A*xl;
      gtov_avg = A*beta0[ii]/Gamma0[ii] + B;
      Vshell = 4.0*M_PI/3.0 * (pow(redge[ii],3) - pow(redge[ii-1],3));
      vedge = (ii < n1) ? 0.5*(vavg[ii-1]+vavg[ii]) : vavg[ii-1];
      gtov[ii] = sqrt(1.0 + pow(vedge/CLIGHT,2) - 2.0*GNEWT*(mr[ii-1] + Vshell*((rhoavg[ii-1]+uavg[ii-1]/SQR(CLIGHT) + eradavg[ii-1])*gtov_avg + vdF[ii-1]))/(redge[ii]*SQR(CLIGHT)));
      A = (gtov[ii] - gtov[ii-1])/dx[0];
      B = gtov[ii-1] - A*xl;
      gtov_avg = A*beta0[ii]/Gamma0[ii] + B;
      mr[ii] = mr[ii-1] + Vshell*((rhoavg[ii-1]+uavg[ii-1]/SQR(CLIGHT) + eradavg[ii-1])*gtov_avg + vdF[ii-1]);
      mstar += Vshell*rhoavg[ii-1];
    }
    if (iter > 0) {
      if (fabs(mlast - mr[n1])/mr[n1] < 1.0e-10) {
        break;
      }
    }
    mlast = mr[n1];
    iter++;
    if (iter > 20) {
      if (mpi_io_proc()) fprintf(stderr,"iter = %d\n", iter);
      fprintf(stderr,"%d %g\n", myrank, mlast);
      fflush(stderr);
      #if (USE_MPI==TRUE)
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
      exit(3);
    }
  }
  if (mpi_io_proc() && iter > 3) fprintf(stderr,"gr monopole required %d iterations\n", iter);

  for (ii=1; ii<n1; ii++) {
    pgasedge = 0.5*(pgasavg[ii-1] + pgasavg[ii]);
    pradedge = 0.5*(pradavg[ii-1] + pradavg[ii]);
    rhoedge  = 0.5*( rhoavg[ii-1] +  rhoavg[ii]);
    uedge    = 0.5*(   uavg[ii-1] +    uavg[ii]);
    grav[ii] = -GNEWT*(mr[ii] + 4.0*M_PI*pow(redge[ii],3)*(pgasedge+pradedge)/SQR(CLIGHT))/pow(redge[ii]*gtov[ii],2) * (rhoedge + (uedge+pgasedge)/SQR(CLIGHT))/rhoedge * dr_dx(startx[0]+ii*dx[0]);
  }
  grav[n1] = -GNEWT*(mr[n1] + 4.0*M_PI*pow(redge[n1],3)*(pgasavg[n1-1]+pradavg[n1-1])/SQR(CLIGHT))/pow(redge[n1]*gtov[n1],2) * (rhoavg[n1-1] + (uavg[n1-1]+pgasavg[n1-1])/SQR(CLIGHT))/rhoavg[n1-1] * dr_dx(startx[0]+n1*dx[0]);

  // now calculate the "GR" potential
  gr_lapse_edge[n1] = -GNEWT*mstar/redge[n1];
  for (ii=n1-1; ii>=0; ii--) {
    gr_lapse_edge[ii] = gr_lapse_edge[ii+1] + 0.5*(grav[ii]+grav[ii+1])*dx[0];   // + because g=-dphi/dr & I'm integrating inwards
  }

  // AS:  ERROR!  This is not threadable, because the new value of grav[ii] depends on the *old* value of grav[ii+1]
  for (ii=0; ii<n1; ii++) {
    xl = startx[0] + ii*dx[0];
    A  = (grav[ii+1]-grav[ii])/dx[0];
    B  = grav[ii] - A*xl;
    dgrav[ii]   = A;
    grav[ii]    = A*beta0[ii]/Gamma0[ii] + B;
    gr_grav[ii] = grav[ii];
  }

  for (ii=0; ii<n1; ii++) {
    xl = startx[0] + ii*dx[0];
    A  = (gr_lapse_edge[ii+1]-gr_lapse_edge[ii])/dx[0];
    B  = gr_lapse_edge[ii] - A*xl;
    gr_lapse[ii] = A*beta0[ii]/Gamma0[ii] + B;
  }

  for (ii=0; ii<n1; ii++) {
    gr_lapse[ii]      = exp(gr_lapse[ii]/SQR(CLIGHT));
    gr_lapse_edge[ii] = exp(gr_lapse_edge[ii]/SQR(CLIGHT));
  }
  gr_lapse_edge[n1] = exp(gr_lapse_edge[n1]/SQR(CLIGHT));

  #else  /* NOT GR_MONOPOLE */
  mpi_global_ireduce_finish(ms, n1, 0);
  mr[0]   = 0.0;
  grav[0] = 0.0;
  for (ii=1; ii<=n1; ii++) {
    mr[ii]   = mr[ii-1] + ms[ii-1];
    grav[ii] = -GNEWT*mr[ii]/SQR(redge[ii])*dr_dx(startx[0]+ii*dx[0]);
  }
  total_mass = mr[n1];
  // AS:  ERROR!  This is not threadable, because the new value of grav[ii] depends on the *old* value of grav[ii+1]
  for (ii=0; ii<n1; ii++) {
#ifdef THIRD_ORDER_MONOPOLE
    xr = startx[0] + (ii+1)*dx[0];
    A  = (grav[ii+2]-grav[ii])/(2.0*dx[0]);
    C  = 0.5*(grav[ii+2] - 2.0*grav[ii+1] + grav[ii])/SQR(dx[0]);
    grav[ii]  = C*alpha0[ii]/Gamma0[ii] + (A-2*C*xr)*beta0[ii]/Gamma0[ii] + (grav[ii+1] - A*xr + C*SQR(xr));
#else
    xl = startx[0] + ii*dx[0];
    A  = (grav[ii+1]-grav[ii])/dx[0];
    B  = grav[ii] - A*xl;
    dgrav[ii] = A;
    grav[ii]  = A*beta0[ii]/Gamma0[ii] + B;
#endif
  }
  #endif  /* ifdef GR_MONOPOLE */

  //ZLOOP {
  GPU_PRAGMA(omp target update to(gr_lapse[:n1],gr_lapse_edge[:n1+1]))
  GPU_PRAGMA(omp target data map(to:grav[:n1+1])) {
  GPU_PRAGMA(omp target teams distribute parallel for)
  for(int II=0;II<cell_count;II++) {
    GET_IJK_FROM_I(II,ii,jj,kk);
    #if (GR_MONOPOLE==TRUE)
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[1] = gr_lapse_edge[ii];
    #endif

    njp  = DJS(ii)/DJS(ii+1);
    jp   = jj*njp;
    vr   = 0.0;
    area = 0.0;
    for (j=jp; j<jp+njp; j++) {
      nkp = DKS(ii,jj)/DKS(ii+1,j);
      kp  = kk*nkp;
      for (k=kp; k<kp+nkp; k++) {
        vr   += NDP_ELEM_LINEAR_F(sim_fdir0,0,ii+1,j,k,0);
        area += ND_ELEM_LINEAR(geom,ii+1,j,k).area[0] + 1.0e-16;
      }
    }
    vr /= area;

    momsrc = NDP_ELEM_LINEAR(p,ii,jj,kk,RHO)*grav[ii];
    NDP_ELEM_LINEAR(sim_src,ii,jj,kk,U1) += momsrc;

    vl = NDP_ELEM_LINEAR_F(sim_fdir0,0,ii,jj,kk,0)/(ND_ELEM_LINEAR(geom,ii,jj,kk).area[0] + 1.0e-16);
    xl = startx[0] + ii*dx[0];
    A  = (vr - vl)/dx[0];
    B  = vl - A*xl;

    NDP_ELEM_LINEAR(sim_src,ii,jj,kk,ETOT) += (A*beta0[ii]/Gamma0[ii] + B)*grav[ii];

    #if (GR_MONOPOLE==TRUE && DO_RADIATION==TRUE)
    GLOOP {
      NDP_ELEM_LINEAR(sim_src,ii,jj,kk,ifrad1+g*NDIM) += gr_lapse[ii]*NDP_ELEM_LINEAR(p,ii,jj,kk,irad1+g)*grav[ii];
    }
    #endif
  }  // end ZLOOP

  // to do: merge into one omp target region.
#if (NDIM==1)
  GPU_PRAGMA(omp target) {
#endif
  #if (GR_MONOPOLE==TRUE)
  ii = istop[0];
  #if (NDIM>1)
  GPU_PRAGMA(omp target teams distribute parallel for)
  JSLOOP(ii,jj) {
  #endif
    #if (NDIM==3)
    KSLOOP(ii,jj,kk) {
    #endif
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[1] = gr_lapse_edge[ii];
    #if (NDIM==3)
    }
    #endif
  #if (NDIM>1)
  }
  #endif

  #if (NDIM==2)
  GPU_PRAGMA(omp target teams distribute parallel for)
  ISLOOP(ii) {
    jj = JS(ii,istop[1]);
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
  }
  #endif

  #if (NDIM==3)
  GPU_PRAGMA(omp target teams distribute parallel for)
  ISLOOP(ii) {
    jj = JS(ii,istop[1]);
    KSLOOP(ii,jj,kk) {
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    }
    JSLOOP(ii,jj) {
      kk = KS(ii,jj,istop[2]);
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    }
  }
  #endif
  #endif /* GR_MONOPOLE */
  }
#if (NDIM==1)
  }
#endif
  TIMER_STOP;
  return;
}

// helper routines for multipole_gravity below

double factorial_ratio(int n, int m)
{
  // returns n!/m!

  double fact = 1.0;

  while (n>m) {
    fact *= n;
    n--;
  }

  while (m>n) {
    fact /= m;
    m--;
  }

  return fact;
}

void get_ylms(double *ylm, double th, double phi, int lmax, int mmax)
{
  int i,l,m;
  double fact,somx2,pmm,norm_fact,smphi,cmphi,norm;
  double cth = cos(th);

  /* AS:  First, compute the associate Legendre polynomial, P[l,m](th) = P[l,m](x), where
   *   x = cos(th), for the non-negative m's.  The negative m's are related as described below.
   * Identities used:
   *   1) P[m,m](x) = (-1)^m * (2*m-1)!! * (1-x^2)^(m/2)
   *   2) P[m+1,m](x) = x * (2*m+1) * P[m,m](x)
   *   3) (l-m+1)*P[l+1,m](x) = (2*l+1)*x*P[l,m](x) - (l+m)*P[l-1,m](x)
   */
  for (m=0; m<=mmax; m++) {
    somx2 = sqrt((1.0-cth)*(1.0+cth));
    fact  = 1.0;
    pmm   = 1.0;
    for (i=1; i<=m; i++) {
      pmm  = -pmm*fact*somx2;
      fact = fact + 2.0;
    }
    ylm[LM_INDEX(m,m)] = pmm;

    if (m != lmax) {
      ylm[LM_INDEX(m+1,m)] = cth*(2.0*m + 1.0)*pmm;

      for (l=m+2; l<=lmax; l++) {
        ylm[LM_INDEX(l,m)] = (cth*(2.0*l-1.0)*ylm[LM_INDEX(l-1,m)]- (l+m-1)*ylm[LM_INDEX(l-2,m)])/(l-m);
      }
    }
  }

  /* AS:  Second, normalize the real spherical harmonics using:
   *               {{  sqrt(2) * sqrt((2*l+1)/(4*pi) * (l-|m|)!/(l+|m|)!) * P[l,|m|](x) * sin(|m|*phi),  for m < 0
   *   Y[l,m](x) = {{            sqrt((2*l+1)/(4*pi))                     * P[l, 0 ](x)               ,  for m = 0
   *               {{  sqrt(2) * sqrt((2*l+1)/(4*pi) * (l- m )!/(l+ m )!) * P[l, m ](x) * cos( m *phi),  for m > 0
   */
  for (m=0; m<=mmax; m++) {
    smphi = sin(m*phi);
    cmphi = cos(m*phi);

    norm_fact = (m == 0) ? 1 : M_SQRT2;
    for (l=m; l<=lmax; l++) {
      norm = norm_fact*sqrt((2.0*l+1.0)/(4.0*M_PI)*factorial_ratio(l-m,l+m));
      ylm[LM_INDEX(l,m)] *= norm;
      if (m>0) {
        ylm[LM_INDEX(l,-m)] = ylm[LM_INDEX(l,m)]*smphi;
        ylm[LM_INDEX(l,m)] *= cmphi;
      }
    }
  }

  return;
}

void get_plms(double *plm, double th, int lmax, int mmax)
{
  int i,l,m;
  double fact,somx2,pmm,norm_fact;
  double cth = cos(th);

  /* AS:  See the first part of get_ylms() above */
  for (m=0; m<=mmax; m++) {
    somx2 = sqrt((1.0-cth)*(1.0+cth));
    fact  = 1.0;
    pmm   = 1.0;
    for (i=1; i<=m; i++) {
      pmm  = -pmm*fact*somx2;
      fact = fact + 2.0;
    }
    plm[LM_INDEX(m,m)] = pmm;

    if (m != lmax) {
      plm[LM_INDEX(m+1,m)] = cth*(2.0*m + 1.0)*pmm;

      for (l=m+2;l<=lmax;l++) {
        plm[LM_INDEX(l,m)] = (cth*(2.0*l-1.0)*plm[LM_INDEX(l-1,m)] - (l+m-1)*plm[LM_INDEX(l-2,m)])/(l-m);
      }
    }

  }

  return;
}

void get_ylm_ints(double *slm, double thmin, double thmax, double phimin, double phimax, int lmax, int mmax)
{
  // integral of Ylm sin\theta over a zone
  int i,l,m;
  double fact,smox2;
  double plmp[LMAX1*LMAX1];
  double plmm[LMAX1*LMAX1];
  double smphi,cmphi;
  double norm_fact,norm;
  double cmin = cos(thmin);
  double cmax = cos(thmax);

  get_plms(plmm, thmin, lmax, mmax);
  get_plms(plmp, thmax, lmax, mmax);

  for (m=0; m<=mmax; m++) {

    if (m>1) {
      /* AS:  Identities used:
      *   P[m,m](x) = (-1)^m * (2*m-1)!! * (1-x^2)^{m/2}
      *   \int sin^{m+1}(th) dth = -sim^m(th)*cos(th)/(m+1) + m/(m+1)*\int sim^{m-1}(th) dth
      */
      // AS:  Fixed the following line
      // Again, some signs got changed, and l's and m's are mixed...
      // slm[LM_INDEX(m,m)] = -(l*(2*l-3)*(2*l-1)*slm[LM_INDEX(m-2,m-2)] + cmax*plmp[LM_INDEX(m,m)] - cmin*plmm[LM_INDEX(m,m)])/(l+1.);
      slm[LM_INDEX(m,m)] = (m*(2*m-1)*(2*m-3)*slm[LM_INDEX(m-2,m-2)] - cmax*plmp[LM_INDEX(m,m)] + cmin*plmm[LM_INDEX(m,m)])/(m+1.0);
    } else if (m==0) {
      /* AS:  P[0,0](th) = 1 */
      slm[LM_INDEX(0,0)] = cmin - cmax;
    } else {
      /* AS:  P[1,1](th) = -sin(th) */
      // slm[LM_INDEX(1,1)] = -0.5*(cmax*sqrt(1-cmax*cmax) + asin(cmax) - cmin*sqrt(1-cmin*cmin) - asin(cmin));
      // Recall that dx = -sin(th) dth; it looks like a sign got dropped...
      slm[LM_INDEX(1,1)] = 0.5*(cmax*sqrt(1-cmax*cmax) - thmax - cmin*sqrt(1-cmin*cmin) + thmin);
    }

    if (m!=lmax) {
      /* AS:  Identities used:
      *   P[m+1,m](x) = (2*m+1)*x*P[m,m](x)
      *   P[m,m](x) = (-1)^m * (2*m-1)!! * (1-x^2)^{m/2}
      */
      // AS:  Fixed the following line
      // slm[LM_INDEX(m+1,m)] = -(sqrt(1-cmax*cmax)*plmp[LM_INDEX(m+1,m+1)] - sqrt(1-cmin*cmin)*plmm[LM_INDEX(m+1,m+1)])/(m+2.);
      slm[LM_INDEX(m+1,m)] = ((1-cmax*cmax)*plmp[LM_INDEX(m,m)] - (1-cmin*cmin)*plmm[LM_INDEX(m,m)])*(2*m+1)/(m+2.0);

      for (l=m+2; l<=lmax; l++) {
        /* AS:  Identities used:
        *   (l-m+1)*P[l+1,m](x) = (2*l+1)*x*P[l,m](x) - (l+m)*P[l-1,m](x)
        *   (1-x^2)*P[l,m](x) = [(l+1)*(l+m)*P[l-1,m](x) - l*(l-m+1)*P[l+1,m](x)]/(2*l+1)
        */
        slm[LM_INDEX(l,m)] = ((l-2)*(l+m-1)*slm[LM_INDEX(l-2,m)] + (2*l-1)*((1-cmax*cmax)*plmp[LM_INDEX(l-1,m)] - (1-cmin*cmin)*plmm[LM_INDEX(l-1,m)]))/((l+1)*(l-m));
      }
    }
  }

  /* AS:  Second, normalize the real spherical harmonics as described in the second part of get_ylms() above.
  *   This will require cmphi = \int sin(m*ph) dph (for m<0) and smphi = \int cos(m*ph) dph (for m>0).
  */
  for (m=0; m<=mmax; m++) {

    if (m==0) {
      smphi = phimax-phimin;
      cmphi = 0.0;
    } else {
      smphi = (sin(m*phimax) - sin(m*phimin))/m;
      cmphi = (cos(m*phimin) - cos(m*phimax))/m;
    }

    norm_fact = (m == 0) ? 1.0 : M_SQRT2;
    for (l=m; l<=lmax; l++) {
      norm = norm_fact*sqrt((2.0*l+1.0)/(4.0*M_PI)*factorial_ratio(l-m,l+m));
      slm[LM_INDEX(l,m)] *= norm;
      if (m>0) {
        slm[LM_INDEX(l,-m)] = slm[LM_INDEX(l,m)]*cmphi;
        slm[LM_INDEX(l,m)] *= smphi;
      } else {
        slm[LM_INDEX(l,0)] *= smphi;
      }
    }
  }

  return;
}

static double *Clm;
static double *Dlm;
static double *ylm;
// static double ND_PTR Phi;
static int mmax, lmax;
#if (USE_LINEAR_ALLOCATION==TRUE)
void multipole_gravity_start(double ** p)
#else
void multipole_gravity_start(double NDP_PTR p)
#endif
{
  static int firstc = 1;
  static int oldstep = 0;
  static int cnt = 0;

  int ii,jj,kk,dd,l,m;
  double rmin,rmax,thmin,thmax,phimin,phimax;
  double rc,rd;
  double x,y,t;

  #if (NDIM==1)
    return;
  #endif

  TIMER_START("multipole_gravity_start");

//  if (multipole_lmax == 0) return;

  if (firstc) {
    ylm = malloc_rank1(   LMAX1*LMAX1, sizeof *ylm);
    Clm = malloc_rank1(n1*LMAX1*LMAX1, sizeof *Clm);
    Dlm = malloc_rank1(n1*LMAX1*LMAX1, sizeof *Dlm);
    // Phi = dendritic_malloc_double();
    firstc = 0;
  }

  for (ii=0; ii<n1*LMAX1*LMAX1; ii++) {
    Clm[ii] = 0.0;
    Dlm[ii] = 0.0;
  }

  #if (NDIM==1)
  lmax = mmax = 0;
  #endif
  #if (NDIM==2)
  lmax = multipole_lmax;
  mmax = 0;
  #endif
  #if (NDIM==3)
  lmax = mmax = multipole_lmax;
  #endif

  ZLOOP {
    rmin = r_of_x(rx_info,startx[0] +  ii   *dx[0]);
    rmax = r_of_x(rx_info,startx[0] + (ii+1)*dx[0]);
    #if (NDIM>1)
    thmin = th_of_x(thx_info,startx[1] +  jj   *DJS(ii)*dx[1]);
    thmax = th_of_x(thx_info,startx[1] + (jj+1)*DJS(ii)*dx[1]);
    #else
    thmin = 0.0;
    thmax = M_PI;
    #endif
    #if (NDIM==3)
    phimin = startx[2] +  kk   *DKS(ii,jj)*dx[2];
    phimax = startx[2] + (kk+1)*DKS(ii,jj)*dx[2];
    #else
    phimin = 0.0;
    phimax = 2.0*M_PI;
    #endif
    get_ylm_ints(ylm, thmin, thmax, phimin, phimax, lmax, mmax);

    for (l=0; l<=lmax; l++) {
      // AS:  The following line is choppy
      //rc = (pow(rmax,l+3) - pow(rmin,l+3))/(l+3);
      rc = pow(rmax,l+3)*(1.0 - pow(rmin/rmax,l+3))/(l+3);
      if      (l==0) rd = 0.5*(SQR(rmax) - SQR(rmin));
      else if (l==1) rd = rmax - rmin;
      else if (l==2) rd = (rmin > 1.0e-16) ? log(rmax/rmin) : 0.0;
      else           rd = (rmin > 1.0e-16) ? (pow(rmin,2-l) - pow(rmax,2-l))/(l-2) : 0.0;
      for (m=-MIN(l,mmax); m<=MIN(l,mmax); m++) {
        Clm[ILM_INDEX(ii,l,m)] += NDP_ELEM_LINEAR(p,ii,jj,kk,RHO)*ylm[LM_INDEX(l,m)]*rc;
        Dlm[ILM_INDEX(ii,l,m)] += NDP_ELEM_LINEAR(p,ii,jj,kk,RHO)*ylm[LM_INDEX(l,m)]*rd;
      }
    }
  }

#if (GR_MONOPOLE==TRUE)
  mpi_global_ireduce_start(Clm, n1*LMAX1*LMAX1,7);
  mpi_global_ireduce_start(Dlm, n1*LMAX1*LMAX1,8);
#else
  mpi_global_ireduce_start(Clm, n1*LMAX1*LMAX1,1);
  mpi_global_ireduce_start(Dlm, n1*LMAX1*LMAX1,2);
#endif

  TIMER_STOP;
  return;
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void multipole_gravity_finish(double ** p)
#else
void multipole_gravity_finish(double NDP_PTR p)
#endif
{
  int ii,jj,kk=0,dd,l,m;
  double Clfac,Dlfac;
  double rm,rp,rc;
  double thm,thp,thc;
  double phm,php,phc;
  double rho,vx,vy,vz,dPhi;
  double vx0,vx1,vy0,vy1,vz0,vz1;
  double dx0,dx1,dx2,dA0,dA1;
  double wx[2],wy[2],wz[2],sgn[2];
  int iii,jjj,kkk,ip,jp,kp,njp,nkp;
  int lev,s;

  #if (NDIM==1)
  return;
  #endif

  TIMER_START("multipole_gravity_finish");

//  if (multipole_lmax == 0) return;

  #if (GR_MONOPOLE==TRUE)
  mpi_global_ireduce_finish(Clm, n1*LMAX1*LMAX1, 7);
  mpi_global_ireduce_finish(Dlm, n1*LMAX1*LMAX1, 8);
  #else
  mpi_global_ireduce_finish(Clm, n1*LMAX1*LMAX1, 1);
  mpi_global_ireduce_finish(Dlm, n1*LMAX1*LMAX1, 2);
  #endif

  /* Compute cumulative sums for C,D coefficients */
  for (ii=1; ii<n1; ii++) {
    for (l=0; l<=lmax; l++) {
      for (m=-MIN(l,mmax); m<=MIN(l,mmax); m++) {
        Clm[ILM_INDEX(ii,l,m)] += Clm[ILM_INDEX(ii-1,l,m)];
      }
    }
  }
  for (ii=n1-2; ii>=0; ii--) {
    for (l=0; l<=lmax; l++) {
      for (m=-MIN(l,mmax); m<=MIN(l,mmax); m++) {
        Dlm[ILM_INDEX(ii,l,m)] += Dlm[ILM_INDEX(ii+1,l,m)];
      }
    }
  }

  /* Compute multipole potential Phi at zone corners */
  thm = 0.5*M_PI;
  phm = M_PI;
  for (ii=istart[0]; ii<=istop[0]; ii++) {
    rm  = r_of_x(rx_info,startx[0]+ ii*dx[0]);

    #if (NDIM>1)
    for (jj=JS(ii,istart[1]); jj<=JS(ii,istop[1]); jj++) {
      thm = th_of_x(thx_info,startx[1] + jj*DJS(ii)*dx[1]);
    #endif

      #if (NDIM==3)
      for (kk=KS(ii,jj,istart[2]); kk<=KS(ii,jj,istop[2]); kk++) {
        phm = startx[2] + kk*DKS(ii,jj)*dx[2];
      #endif

        get_ylms(ylm, thm, phm, lmax, mmax);

        ND_ELEM_LINEAR(sim_Phi,ii,jj,kk) = 0.0;

        if (rm < 1.0e-16) {
        } else if (ii==n1) {
          for (l=1; l<=lmax; l++) {
            Clfac = -4.0*M_PI*GNEWT/(2*l+1)/pow(rm,l+1);
            for (m=-MIN(l,mmax); m<=MIN(l,mmax); m++) {
              ND_ELEM_LINEAR(sim_Phi,ii,jj,kk) += Clfac*ylm[LM_INDEX(l,m)]*Clm[ILM_INDEX(ii-1,l,m)];
            }
          }
        } else {
          for (l=1; l<=lmax; l++) {
            Clfac = -4.0*M_PI*GNEWT/(2*l+1)/pow(rm,l+1);
            Dlfac = -4.0*M_PI*GNEWT/(2*l+1)*pow(rm,l  );
            for (m=-MIN(l,mmax); m<=MIN(l,mmax); m++) {
              ND_ELEM_LINEAR(sim_Phi,ii,jj,kk) += Clfac*ylm[LM_INDEX(l,m)]*Clm[ILM_INDEX(ii-1,l,m)];
              ND_ELEM_LINEAR(sim_Phi,ii,jj,kk) += Dlfac*ylm[LM_INDEX(l,m)]*Dlm[ILM_INDEX(ii  ,l,m)];
            }
          }
        }

      #if (NDIM==3)
      }
      #endif
    #if (NDIM>1)
    }
    #endif
  }

  /* Compute gravity source terms:  For the component of grad Phi in a given direction, need to loop
   *   over each corner of a given zone and compute grad Phi along edges in that direction, then
   *   perform a weighted sum of grad Phi along each edge to get estimate the volume average.  This
   *   is a bit simpler if we use weights, where the weight for each lower edge uses the fractional
   *   gradient at the upper part of the zone and vice versa.  Note that the weights in a given
   *   direction always add to 1.  The sgn[] weight, which is either +1 for outer or -1 for inner, is
   *   used to perform the differencing inside the loops.  Finally, <rho*v> is computed using the mass
   *   fluxes and dividing by the area, then averaging using the same weights as before.
   */
  sgn[0] = -1.0;  sgn[1] = 1.0;
  wy[0] = 0.5;  wy[1] = 0.5;
  wz[0] = 0.5;  wz[1] = 0.5;
  dx0 = 1.0;  dx1 = 1.0;  dx2 = 1.0;

  ISLOOP(ii) {
    rm  = r_of_x(rx_info,startx[0]+ ii   *dx[0]);
    rp  = r_of_x(rx_info,startx[0]+(ii+1)*dx[0]);
    /* Use beta/Gamma to get <x[0]> */
    rc  = r_of_x(rx_info,beta0[ii]/Gamma0[ii]);
    wx[0] = (rp-rc)/(rp-rm);
    wx[1] = (rc-rm)/(rp-rm);
    dx0 = dx[0];

    /* Compute the dendritic level in jj for this ii */
  	lev = 0;
  	s   = DJS(ii);
  	while (s >>= 1) lev++;

    #if (NDIM>1)
    JSLOOP(ii,jj) {
      thm = th_of_x(thx_info,startx[1] +  jj   *DJS(ii)*dx[1]);
      thp = th_of_x(thx_info,startx[1] + (jj+1)*DJS(ii)*dx[1]);
      /* Use beta/Gamma to get <x[1]> */
      thc = th_of_x(thx_info,beta1s[lev][jj]/Gamma1s[lev][jj]);
      wy[0] = (thp-thc)/(thp-thm);
      wy[1] = (thc-thm)/(thp-thm);
      dx1 = DJS(ii)*dx[1];
    #endif
      #if (NDIM==3)
      KSLOOP(ii,jj,kk) {
        phm = startx[2] +  kk   *DKS(ii,jj)*dx[2];
        php = startx[2] + (kk+1)*DKS(ii,jj)*dx[2];
        /* <x[2]> is just the arithemetic mean in the 2-direction */
        phc = 0.5*(phm + php);
        wz[0] = (php-phc)/(php-phm);
        wz[1] = (phc-phm)/(php-phm);
        dx2 = DKS(ii,jj)*dx[2];
      #endif

        rho = NDP_ELEM_LINEAR(p,ii,jj,kk,RHO);

        /* Compute grad Phi in 0-direction */
        dPhi = 0.0;
        for (iii=0; iii<=1; iii++) {
          ip = ii+iii;
          #if (NDIM>1)
          for (jjj=0; jjj<=1; jjj++) {
            jp = (jj+jjj)*DJS(ii)/DJS(ip);
          #endif
            #if (NDIM==3)
            for (kkk=0; kkk<=1; kkk++) {
              kp = (kk+kkk)*DKS(ii,jj)/DKS(ip,jp);
            #endif
              dPhi  += ND_ELEM_LINEAR(sim_Phi,ip,jp,kp)*sgn[iii]*wy[jjj]*wz[kkk];
            #if (NDIM==3)
            }
            #endif
          #if (NDIM>1)
          }
          #endif
        }

        /* Compute rho*v in the 0-direction using mass fluxes */
        dA0 = ND_ELEM_LINEAR(geom,ii,jj,kk).area[0] + 1.0e-16;
        vx0 = NDP_ELEM_LINEAR_F(sim_fdir0,0,ii,jj,kk,RHO)/dA0;
        njp = DJS(ii)/DJS(ii+1);
        jp = jj*njp;
        dA1 = 0.0;
        vx1 = 0.0;
        for (jjj=jp; jjj<jp+njp; jjj++) {
          nkp = DKS(ii,jj)/DKS(ii+1,jjj);
          kp = kk*nkp;
          for (kkk=kp; kkk<kp+nkp; kkk++) {
            dA1 += ND_ELEM_LINEAR(geom,ii+1,jjj,kkk).area[0] + 1.0e-16;
            vx1 += NDP_ELEM_LINEAR_F(sim_fdir0,0,ii+1,jjj,kkk,RHO);
          }
        }
        vx1 /= dA1;
        vx = vx0*wx[0] + vx1*wx[1];

        NDP_ELEM_LINEAR(sim_src,ii,jj,kk,U1  ) -= rho*dPhi/dx0;
        NDP_ELEM_LINEAR(sim_src,ii,jj,kk,ETOT) -= vx*dPhi/dx0;

        #if (NDIM>1)
        /* Compute grad Phi in 1-direction */
        dPhi = 0.0;
        for (iii=0; iii<=1; iii++) {
          ip = ii+iii;
          #if (NDIM>1)
          for (jjj=0; jjj<=1; jjj++) {
            jp = (jj+jjj)*DJS(ii)/DJS(ip);
          #endif
            #if (NDIM==3)
            for (kkk=0; kkk<=1; kkk++) {
              kp = (kk+kkk)*DKS(ii,jj)/DKS(ip,jp);
            #endif
              dPhi  += ND_ELEM_LINEAR(sim_Phi,ip,jp,kp)*wx[iii]*sgn[jjj]*wz[kkk];
            #if (NDIM==3)
            }
            #endif
          #if (NDIM>1)
          }
          #endif
        }

        /* Compute rho*v in the 1-direction using mass fluxes (recall that the flux of rho is rho*v*dA) */
        /* Outer 1-edge */
        dA1 = 0.0;
        vy1 = 0.0;
        if (DKS(ii,jj+1) > DKS(ii,jj)) {
          // Coarsening boundary, 1 outer neighbor, shifted area index, half area, flux index equals original index
          dA1 += 0.5*ND_ELEM_LINEAR(geom,ii,jj+1,kk/2).area[1] + 1.0e-16;
          vy1 += NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,kk,RHO);
        } else if (DKS(ii,jj+1) < DKS(ii,jj)) {
          // Refinement boundary, 2 outer neighbors, shifted area index, original area, flux index equals area index
          dA1 += ND_ELEM_LINEAR(geom,ii,jj+1,2*kk).area[1] + ND_ELEM_LINEAR(geom,ii,jj+1,2*kk+1).area[1] + 1.0e-16;
          vy1 += NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,2*kk,RHO) + NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,2*kk+1,RHO);
        } else {
          // Regular boundary, 1 outer neighbor, original area index, original area, flux index equals original index
          dA1 += ND_ELEM_LINEAR(geom,ii,jj+1,kk).area[1] + 1.0e-16;
          vy1 += NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,kk,RHO);
        }
        vy1 /= dA1;

        /* Inner 1-edge */
        dA0 = 0.0;
        vy0 = 0.0;
        if (DKS(ii,jj-1) < DKS(ii,jj)) {
          // Coarsening boundary, 2 inner neighbors, original area index, half area, shifted flux index
          dA0 += ND_ELEM_LINEAR(geom,ii,jj,kk).area[1] + 1.0e-16;
          vy0 += NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj,2*kk,RHO) + NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj,2*kk+1,RHO);
        } else {
          // Refinement boundary, 1 inner neighbor, original area index, original area, flux index equals area index
          // Regular boundary, 1 inner neighbor, original area index, original area, flux index equals original index
          dA0 += ND_ELEM_LINEAR(geom,ii,jj,kk).area[1] + 1.0e-16;
          vy0 += NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj,kk,RHO);
        }
        vy0 /= dA0;
        vy = vy0*wy[0] + vy1*wy[1];


        NDP_ELEM_LINEAR(sim_src,ii,jj,kk,U2  ) -= rho*dPhi/dx1;
        NDP_ELEM_LINEAR(sim_src,ii,jj,kk,ETOT) -= vy*dPhi/dx1;
        #endif

        #if (NDIM==3)
        /* Compute grad Phi in 2-direction */
        dPhi = 0.0;
        for (iii=0; iii<=1; iii++) {
          ip = ii+iii;
          #if (NDIM>1)
          for (jjj=0; jjj<=1; jjj++) {
            jp = (jj+jjj)*DJS(ii)/DJS(ip);
          #endif
            #if (NDIM==3)
            for (kkk=0; kkk<=1; kkk++) {
              kp = (kk+kkk)*DKS(ii,jj)/DKS(ip,jp);
            #endif
              dPhi  += ND_ELEM_LINEAR(sim_Phi,ip,jp,kp)*wx[iii]*wy[jjj]*sgn[kkk];
            #if (NDIM==3)
            }
            #endif
          #if (NDIM>1)
          }
          #endif
        }

        /* Compute rho*v in the 2-direction using mass fluxes */
        dA0 = ND_ELEM_LINEAR(geom,ii,jj,kk  ).area[2] + 1.0e-16;
        vz0 = NDP_ELEM_LINEAR_F(sim_fdir2,2,ii,jj,kk  ,RHO)/dA0;
        dA1 = ND_ELEM_LINEAR(geom,ii,jj,kk+1).area[2] + 1.0e-16;
        vz1 = NDP_ELEM_LINEAR_F(sim_fdir2,2,ii,jj,kk+1,RHO)/dA1;
        vz  = vz0*wz[0] + vz1*wz[1];

        NDP_ELEM_LINEAR(sim_src,ii,jj,kk,U3  ) -= rho*dPhi/dx2;
        NDP_ELEM_LINEAR(sim_src,ii,jj,kk,ETOT) -= vz*dPhi/dx2;
        #endif

      #if (NDIM==3)
      }
      #endif
    #if (NDIM>1)
    }
    #endif
  }

  TIMER_STOP;
  return;
}

#endif

#if (USE_LINEAR_ALLOCATION==TRUE)
void gravity_start(double ** p)
#else
void gravity_start(double NDP_PTR p)
#endif
{
  TIMER_START("gravity_start");

  #if (GRAV==USER_GRAV)
  user_gravity(p);
  #endif

  #if (GRAV==FIXED_GRAV)
  fixed_gravity(p);
  #endif

  #if (GRAV==SPHERICAL_MONOPOLE_GRAV)
  monopole_gravity_start(p);
  #endif

  #if (GRAV==SPHERICAL_MULTIPOLE_GRAV)
  monopole_gravity_start(p);
  multipole_gravity_start(p);
  #endif

  TIMER_STOP;
  return;
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void gravity_finish(double ** p)
#else
void gravity_finish(double NDP_PTR p)
#endif
{
  TIMER_START("gravity_finish");

  #if (GRAV==PRESCRIBED_GRAV)
  prescribed_gravity(p);
  #endif

  #if (GRAV==SPHERICAL_MONOPOLE_GRAV)
  monopole_gravity_finish(p);
  #endif

  #if (GRAV==SPHERICAL_MULTIPOLE_GRAV)
  monopole_gravity_finish(p);
  multipole_gravity_finish(p);
  #endif

  TIMER_STOP;
  return;
}
