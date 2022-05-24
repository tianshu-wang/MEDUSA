#include "../decs.h"

#if ((THCOORD==MU) || (THCOORD==NICEMU))
#error Theta coord options NICEMU and MU are broken
#endif

static double r_;
#if(RCOORD==SINH_MODIFIED)
//static double rx_info.facta=10*1e5,rx_info.factb=15*1e5,rx_info.factc=5*1e5;
#endif
double muBnorm;
double muB = 10.1;

// double **mu1s,**nu1s,**alpha1s,**beta1s,**Gamma1s,**alpha2s,**beta2s,**Gamma2s;

double rtrans_solve(double dx, double rout)
{
  double rs,fl,fr,fm,fp,f;
  double rg = 100*rout;


#if(RCOORD!=UNIFORM)
  // do a little bisection first
  rs = 0.1*dx;
  rx_info.rtrans = rs;
  do {
    rs *= 1.1;
    rx_info.rtrans = rs;
    fl = r_of_x(rx_info,n1*dx)-rout;
  } while(isinf(fl) || isnan(fl));
  rx_info.rtrans = rg;
  fr = r_of_x(rx_info,n1*dx)-rout;
  if (fl*fr > 0.0) {
    fprintf(stderr,"did not bracket rtrans in rtrans_solve %g->%g  %g->%g\n", rs, fl, rg, fr);
    exit(1);
  }

  #if(RCOORD==SINH)
  do {
    rx_info.rtrans = 0.5*(rs + rg);
    fm = r_of_x(rx_info,n1*dx)-rout;
    if (fm*fl <= 0) {
      fr = fm;
      rg = rx_info.rtrans;
    } else {
      fl = fm;
      rs = rx_info.rtrans;
    }
  } while(fabs((rs-rg)/rg) > 1e-10);
  //} while(fabs((rs-rg)/rg) > 0.1);
  rg = 0.5*(rs+rg);
//  do {
//    rs = rg;
//    rx_info.rtrans = rs;
//    f = r_of_x(rx_info,n1*dx) - rout;
//    fp = sinh(n1*dx/rg) - n1*dx*cosh(n1*dx/rg)/rg;
//    rg -= f/fp;
//  } while(fabs((rs-rg)/rg) > 1.0e-10);
  #endif
  #if(RCOORD==SINH_MODIFIED)
  do {
    rx_info.rtrans = 0.5*(rs + rg);
    fm = r_of_x(rx_info,n1*dx)-rout;
    if (fm*fl <= 0) {
      fr = fm;
      rg = rx_info.rtrans;
    } else {
      fl = fm;
      rs = rx_info.rtrans;
    }
  } while(fabs((rs-rg)/rg) > 1e-10);
  rg = 0.5*(rs+rg);
  #endif
#endif
  return fabs(rg);
}

void ijk_to_x(int i, int j, int k, double x[SPACEDIM])
{
  x[0] = startx[0] + (i+0.5)*dx[0];
  #if (NDIM>1)
  x[1] = startx[1] + (j+0.5)*DJS(i)*dx[1];
  #else
  x[1] = 0.5;
  #endif
  #if (NDIM==3)
  x[2] = startx[2] + (k+0.5)*DKS(i,j)*dx[2];
  #else
  x[2] = 0.0;
  #endif

  return;
}

void x_to_ijk(double x[], int index[])
{
  index[0] = (x[0]-startx[0])/dx[0];
  #if(NDIM>1)
  index[1] = (x[1]-startx[1])/(DJS(index[0])*dx[1]);
  #endif
  #if(NDIM==3)
  index[2] = (my_mod(x[2],n3*dx[2])-startx[2])/(DKS(index[0],index[1])*dx[2]);
  #endif
}

void ijk_to_Cart(int i, int j, int k, double *x)
{
  double r   = r_of_x(rx_info,startx[0] + i*dx[0]);
  double th  = th_of_x(thx_info,startx[1] + j*dx[1]);
  double phi = startx[2] + k*dx[2];
  double sth = sin(th);

  #if (NDIM==1)
  x[0] = r;
  #elif (NDIM==2)
  x[0] = r*sth;
  x[1] = r*cos(th);
  #elif (NDIM==3)
  x[0] = r*sth*cos(phi);
  x[1] = r*sth*sin(phi);
  x[2] = r*cos(th);
  #endif

  return;
}

// ************************************* //
//                                       //
//       basic coordinate mappings       //
//                                       //
// ************************************* //
double eps_of_x(double x) {
  return r_of_x(rx_info,x)/r_ - 1.0;
}

double deps_dx(double x) {
  return dr_dx(x)/r_;
}

double comp_x_from_r(double r, double xlow, double xhigh) {
  r_ = r;
  return dynamic_root_find(eps_of_x, deps_dx, 0.5*(xlow+xhigh), xlow, xhigh, 1e-10, 1000);
}

double sigmoid(double x){
  return 1.0/(1.0+exp(-x));
}

double r_of_x(r_of_x_struct_t rx_info,double x)
{
  #if (RCOORD==UNIFORM)
  return fabs(x);
  #endif

  #if (RCOORD==SINH)
  return rx_info.rtrans*sinh(fabs(x)/rx_info.rtrans);
  #endif

  #if (RCOORD==SINH_MODIFIED)
  return rx_info.rtrans*sinh(fabs(x)/rx_info.rtrans)-rx_info.facta*sigmoid((fabs(x)-rx_info.factb)/rx_info.factc)+rx_info.facta*sigmoid(-rx_info.factb/rx_info.factc);
  #endif

  #if (RCOORD==EXPSINH)
  return 2.0*(rsparse_fact-1.0)*rsparse_trans*(1.0 - exp(-fabs(x)/rsparse_trans)) + rx_info.rtrans*sinh(fabs(x)/rx_info.rtrans);
  #endif
}

double x_of_r(double r)
{
  #if (RCOORD==UNIFORM)
  return r;
  #elif (RCOORD==SINH)
  double y = r/rx_info.rtrans;
  return rx_info.rtrans*log(y + sqrt(SQR(y) + 1.0));
  #else
  return comp_x_from_r(r, 0.0, n1*ccsn_dr_min);
  #endif
}

double dr_dx(double x)
{
  #if (RCOORD==UNIFORM)
  return 1.0;
  #endif

  #if (RCOORD==SINH)
  return cosh(fabs(x)/rx_info.rtrans);
  #endif

  #if (RCOORD==SINH_MODIFIED)
  double expval = exp(-(fabs(x)-rx_info.factb)/rx_info.factc);
  return cosh(fabs(x)/rx_info.rtrans)-rx_info.facta*expval/rx_info.factc/(1+expval)/(1+expval);
  #endif

  #if (RCOORD==EXPSINH)
  return 2.0*(rsparse_fact-1.0)*exp(-fabs(x)/rsparse_trans) + cosh(fabs(x)/rx_info.rtrans);
  #endif
}

double d2r_dx2(double x)
{
  #if (RCOORD==UNIFORM)
  return 0.0;
  #endif

  #if (RCOORD==SINH)
  return sinh(fabs(x)/rx_info.rtrans)/rx_info.rtrans;
  #endif

  #if (RCOORD==SINH_MODIFIED)
  double expval = exp(-(fabs(x)-rx_info.factb)/rx_info.factc);
  return sinh(fabs(x)/rx_info.rtrans)/rx_info.rtrans-rx_info.facta/rx_info.factc/rx_info.factc*(2*expval*expval/(1+expval)/(1+expval)/(1+expval)-expval/(1+expval)/(1+expval));
  #endif

  #if (RCOORD==EXPSINH)
  return -2.0*(rsparse_fact-1.0)/rsparse_trans*exp(-fabs(x)/rsparse_trans) + sinh(fabs(x)/rx_info.rtrans)/rx_info.rtrans;
  #endif
}

// muBnorm = 1./asin(1/sqrt(muB));

double th_of_x(th_of_x_struct_t thx_info,double x)
{
  // Reflect the coordinates about the polar axis
  if (x < -1.0) x = -2.0 - x;
  if (x >  1.0) x =  2.0 - x;

  #if (THCOORD==UNIFORM)
  return 0.5*M_PI*(x+1.0);
  #endif

  #if (THCOORD==POLYTH)
  return thx_info.poly_norm*x*(1.0 + pow(x/thx_info.poly_xt,thx_info.poly_alpha)/(thx_info.poly_alpha+1.0)) + 0.5*M_PI;
  #endif

  #if (THCOORD==NICEMU)
  return thx_info.nice_norm*(asin(1.0/sqrt(thx_info.nice_alpha)) + atan(x/sqrt(thx_info.nice_alpha - SQR(x))));
  #endif

  #if (THCOORD==MU)
  return acos(-x);
  #endif
}

double x_of_th(double th)
{
  #if (THCOORD==UNIFORM)
  return 2.0*th/M_PI - 1.0;
  #elif (THCOORD==NICEMU)
  return sqrt(thx_info.nice_alpha)*tan((-asin(1.0/sqrt(thx_info.nice_alpha))*thx_info.nice_norm + th)/thx_info.nice_norm)/sqrt(1.0 + pow(tan((-asin(1.0/sqrt(thx_info.nice_alpha))*thx_info.nice_norm + th)/thx_info.nice_norm),2.0));
  #elif (THCOORD==MU)
  return -cos(th);
  #else
  return 1.0;
  #endif
}

double dth_dx(double x)
{
  // Reflect the coordinates about the polar axis
  if (x < -1.0) x = -2.0 - x;
  if (x >  1.0) x =  2.0 - x;

  #if (THCOORD==UNIFORM)
  return 0.5*M_PI;
  #endif

  #if (THCOORD==POLYTH)
  return thx_info.poly_norm*(1.0 + pow(x/thx_info.poly_xt,thx_info.poly_alpha));
  #endif

  #if (THCOORD==NICEMU)
  return thx_info.nice_norm/sqrt(thx_info.nice_alpha - SQR(x));
  #endif

  #if (THCOORD==MU)
  return 1.0/sqrt(1.0 - SQR(x) + 1.0e-16);
  #endif
}

double d2th_dx2(double x)
{
  // Reflect the coordinates about the polar axis
  if (x < -1.0) x = -2.0 - x;
  if (x >  1.0) x =  2.0 - x;

  #if (THCOORD==UNIFORM)
  return 0.0;
  #endif

  #if (THCOORD==POLYTH)
  return thx_info.poly_norm*thx_info.poly_alpha*pow(x/thx_info.poly_xt,thx_info.poly_alpha-1.0)/thx_info.poly_xt;
  #endif

  #if (THCOORD==NICEMU)
  return thx_info.nice_norm*x/pow(thx_info.nice_alpha - SQR(x),1.5);
  #endif

  #if (THCOORD==MU)
  return x/pow(1.0-SQR(x),1.5);
  #endif
}


// ************************************* //
//                                       //
//     required coordinate integrals     //
//                                       //
// ************************************* //

// ALL DONE NUMERICALLY FOR CONVENIENCE AND GENERALITY

// radial integrals first

double r2dr_dx(double x)
{
  return pow(r_of_x(rx_info,x),2)*dr_dx(x);
}

double Vr(double x)
{
  // \int_x^x+dx r^2  dr/dx dx
  return romb(x,x+dx[0],r2dr_dx);
}

double r2d2r_dx2(double x)
{
  return pow(r_of_x(rx_info,x),2)*d2r_dx2(x);
}

double int_r2d2r_dx2(double x)
{
  // \int_x^x+dx r^2 d^2r/dx^2 dx
  return romb(x,x+dx[0],r2d2r_dx2);
}

double r3(double x)
{
  return pow(r_of_x(rx_info,x),3);
}

double int_r3(double x)
{
  // \int_x^x+dx r^3 dx
  return romb(x,x+dx[0],r3);
}

double r2dr_dx3(double x)
{
  return pow(r_of_x(rx_info,x),2)*pow(dr_dx(x),3);
}

double int_r2dr_dx3(double x)
{
  // \int_x^x+dx r^2 (dr/dx)^3 dx
  return romb(x,x+dx[0],r2dr_dx3);
}

double rdr_dx2(double x)
{
  return r_of_x(rx_info,x)*pow(dr_dx(x),2);
}

double int_rdr_dx2(double x)
{
  // \int_x^x+dx r (dr/dx)^2 dx
  return romb(x,x+dx[0],rdr_dx2);
}

double xr1g(double x)
{
  return x*pow(r_of_x(rx_info,x),2)*dr_dx(x);
}

double xr2g(double x)
{
  return x*x*pow(r_of_x(rx_info,x),2)*dr_dx(x);
}

double xr3g(double x)
{
  return pow(x,3)*pow(r_of_x(rx_info,x),2)*dr_dx(x);
}

double xr4g(double x)
{
  return pow(x,4)*pow(r_of_x(rx_info,x),2)*dr_dx(x);
}

// now do angular integrals

double sthdth_dx(double x)
{
  return sin(th_of_x(thx_info,x))*dth_dx(x);
}

double sthdth_dx3(double x)
{
  return sin(th_of_x(thx_info,x))*pow(dth_dx(x),3);
}

double sth3dth_dx(double x)
{
  return pow(sin(th_of_x(thx_info,x)),3)*dth_dx(x);
}

double sthd2th_dx2(double x)
{
  return sin(th_of_x(thx_info,x))*d2th_dx2(x);
}

double sth2cth(double x)
{
  double th = th_of_x(thx_info,x);
  return pow(sin(th),2) * cos(th);
}

double cthdth_dx2(double x)
{
  return cos(th_of_x(thx_info,x))*pow(dth_dx(x),2);
}

double xth1g(double x)
{
  return x*sin(th_of_x(thx_info,x))*dth_dx(x);
}

double xth2g(double x)
{
  return x*x*sin(th_of_x(thx_info,x))*dth_dx(x);
}

double xth3g(double x)
{
  return pow(x,3)*sin(th_of_x(thx_info,x))*dth_dx(x);
}

double xth4g(double x)
{
  return pow(x,4)*sin(th_of_x(thx_info,x))*dth_dx(x);
}

// ************************************* //
//                                       //
//     finished coordinate integrals     //
//                                       //
// ************************************* //


void x_to_rthphi(double x[SPACEDIM], double rvec[SPACEDIM])
{
  rvec[0] = r_of_x(rx_info,x[0]);
  rvec[1] = th_of_x(thx_info,x[1]);
  rvec[2] = x[2];
}

void ijk_to_rthphi(int i, int j, int k, double rvec[SPACEDIM])
{
  double x[SPACEDIM];

  ijk_to_x(i, j, k, x);

  rvec[0] = r_of_x(rx_info,x[0]);
  rvec[1] = th_of_x(thx_info,x[1]);
  rvec[2] = x[2];

  return;
}

void ijk_to_r_dr(int i, int j, int k, double r[])
{
  double x[SPACEDIM];

  ijk_to_x(i, j, k, x);

  r[0] = r_of_x(rx_info,x[0]-0.5*dx[0]);
  r[1] = r_of_x(rx_info,x[0]          );
  r[2] = r_of_x(rx_info,x[0]+0.5*dx[0]);
}


double ijk_to_r(int i, int j, int k, double rhat[SPACEDIM])
{
  double r,x[SPACEDIM];

  ijk_to_x(i, j, k, x);
  x[0] = beta0[i]/Gamma0[i];
  r = r_of_x(rx_info,x[0]);
  rhat[0] = dr_dx(x[0]);
  rhat[1] = 0.0;
  rhat[2] = 0.0;

  return r;
}


void vec_transform_to_xyz(double *vp, int i, int j, int k)
{
  // A.S.  This is only ever used in io.c in code that ignores the
  // dendritic grid.  Therefore, it's okay that it uses *dx[1] and *dx[2]
  int l,m;
  double x[SPACEDIM];
  double v[SPACEDIM];
  double lam[SPACEDIM][SPACEDIM];
  double r,th,sth,cth,sph,cph;
  double drdx,dthdx;

  x[0]  = startx[0] + (i+0.5)*dx[0];
  r     = r_of_x(rx_info,x[0]);

  x[1]  = startx[1] + (j+0.5)*dx[1];
  th    = th_of_x(thx_info,x[1]);
  sth   = sin(th);
  cth   = cos(th);

  #if (NDIM==3)
  x[2]  = startx[2] + (k+0.5)*dx[2];
  sph   = sin(x[2]);
  cph   = cos(x[2]);
  #else
  cph   = 1.0;
  sph   = 0.0;
  #endif

  lam[0][0] = sth*cph;
  lam[0][1] = cth*cph;
  lam[0][2] = -sph;

  lam[1][0] = sth*sph;
  lam[1][1] = cth*sph;
  lam[1][2] = cph;

  lam[2][0] = cth;
  lam[2][1] = -sth;
  lam[2][2] = 0.0;

  v[0] = v[1] = v[2] = 0.0;

  for (l=0; l<SPACEDIM; l++) {
    for (m=0; m<SPACEDIM; m++) {
      v[l] += lam[l][m]*vp[m];
    }
  }

  #if (NDIM==2)
  vp[0] = v[0];
  vp[1] = v[2];
  vp[2] = v[1];
  #else
  vp[0] = v[0];
  vp[1] = v[1];
  vp[2] = v[2];
  #endif

  return;
}

void init_coords()
{
  int ii;

  rcenter = malloc_rank1(n1+2*NG  , sizeof *rcenter);
  redge   = malloc_rank1(n1+2*NG+1, sizeof *redge  );
  dr      = malloc_rank1(n1+2*NG  , sizeof *dr     );

  thx_info.nice_norm = M_PI/(asin(1.0/sqrt(thx_info.nice_alpha)) + atan(1.0/sqrt(thx_info.nice_alpha - 1.0)));

  thx_info.poly_norm = 0.5*M_PI*1.0/(1.0 + 1.0/(thx_info.poly_alpha+1.0) * 1.0/pow(thx_info.poly_xt,thx_info.poly_alpha));

  #if (NDIM<3)
  dx[2] = 2.0*M_PI;
  #endif

  #if (NDIM<2)
  dx[1] = 2.0;
  #endif

  rcenter += NG;
  redge   += NG;
  dr      += NG;

  for (ii=-NG; ii<n1+NG; ii++) {
    redge[ii]   = r_of_x(rx_info,fabs( ii   *dx[0] + startx[0]));
    dr[ii]      = r_of_x(rx_info,fabs((ii+1)*dx[0] + startx[0])) - redge[ii];
    rcenter[ii] = redge[ii] + 0.5*dr[ii];
    //if(myrank==0){
    //  fprintf(stderr,"%d %e\n",ii,rcenter[ii]);
    //}
  }
}

void init_interp_vol()
{
  int i,j,k,lev,s,n2s,n3s;
  double x0,x1,r0,r1;

  //#pragma omp parallel for private(i,x0,x1)
  for (i=-NG; i<n1+NG; i++) {
    x0 = startx[0] +  i   *dx[0];
    x1 = startx[0] + (i+1)*dx[0];
    Gamma0[i] = romb(x0,x1,r2dr_dx);
    beta0 [i] = romb(x0,x1,xr1g   );
    alpha0[i] = romb(x0,x1,xr2g   );
  }

  #if (NDIM>1)
  // now deal with dendritic grid stuff
  Gamma1s = malloc_rank1(max_jrefine_level+1, sizeof *Gamma1s);
  beta1s  = malloc_rank1(max_jrefine_level+1, sizeof *beta1s );
  alpha1s = malloc_rank1(max_jrefine_level+1, sizeof *alpha1s);

  for (lev=0; lev<=max_jrefine_level; lev++) {
    //s = 1<<lev;//s=2^lev
    s = 1;
    n2s = global_grid_dims[1]/s + 2*NG;
    Gamma1s[lev]  = malloc_rank1(n2s, sizeof *Gamma1s[lev]);
    Gamma1s[lev] += NG; //- istart[1]/s;
    beta1s [lev]  = malloc_rank1(n2s, sizeof *beta1s [lev]);
    beta1s [lev] += NG; //- istart[1]/s;
    alpha1s[lev]  = malloc_rank1(n2s, sizeof *alpha1s[lev]);
    alpha1s[lev] += NG; //- istart[1]/s;
  }

  //#pragma omp parallel for private(lev,s,j,x0,x1)
  for (lev=0; lev<=max_jrefine_level; lev++) {
    s = 1<<lev;
    for (j=istart[1]/s-NG; j<istop[1]/s+NG; j++) {
      x0 = startx[1] +  j   *s*dx[1];
      x1 = startx[1] + (j+1)*s*dx[1];
      Gamma1s[lev][j] = romb(x0,x1,sthdth_dx);
      beta1s [lev][j] = romb(x0,x1,xth1g    );
      alpha1s[lev][j] = romb(x0,x1,xth2g    );
    }
  }
  #endif

  #if (NDIM>2)
  // now deal with dendritic grid stuff
  Gamma2s = malloc_rank1(max_krefine_level+1, sizeof *Gamma2s);
  beta2s  = malloc_rank1(max_krefine_level+1, sizeof *beta2s );
  alpha2s = malloc_rank1(max_krefine_level+1, sizeof *alpha2s);

  for (lev=0; lev<=max_krefine_level; lev++) {
    //s = 1<<lev;
    s = 1;
    n3s = my_grid_dims[2]/s + 2*NG;
    Gamma2s[lev]  = malloc_rank1(n3s, sizeof *Gamma2s[lev]);
    Gamma2s[lev] += NG; //- istart[2]/s;
    beta2s [lev]  = malloc_rank1(n3s, sizeof *beta2s [lev]);
    beta2s [lev] += NG; //- istart[2]/s;
    alpha2s[lev]  = malloc_rank1(n3s, sizeof *alpha2s[lev]);
    alpha2s[lev] += NG; //- istart[2]/s;
  }

  //#pragma omp parallel for private(lev,s,k,x0,x1)
  for (lev=0; lev<=max_krefine_level; lev++) {
    s = 1<<lev;
    for (k=istart[2]/s-NG; k<istop[2]/s+NG; k++) {
      x0 = startx[2] +  k   *s*dx[2];
      x1 = startx[2] + (k+1)*s*dx[2];
      Gamma2s[lev][k] = s*dx[2];
      beta2s [lev][k] = ( SQR(x1) -  SQR(x0))/2.0;
      alpha2s[lev][k] = (CUBE(x1) - CUBE(x0))/3.0;
    }
  }
  #endif

  return;
}

void calc_interp_coeffs1(int i, double *alpha, double *beta, double *Gamma, double *x)
{
  int j,jstart,lev,s;

  lev = 0;
  s = DJS(i);
  while (s >>= 1) lev++;
  jstart = JS(i,istart[1]);
  JSGLOOP(i,j) {
    x    [j-jstart] = startx[1] + j*DJS(i)*dx[1];
    Gamma[j-jstart] = Gamma1s[lev][j];
    beta [j-jstart] = beta1s [lev][j];
    alpha[j-jstart] = alpha1s[lev][j];
  }

  return;
}

void calc_interp_coeffs2(int i, int j, double *alpha, double *beta, double *Gamma, double *x)
{
  int k,kstart,kstop,lev,s;

  lev = 0;
  s = DKS(i,j);
  while (s >>= 1) lev++;
  kstart = KS(i,j,istart[2]);
  KSGLOOP(i,j,k) {
    x    [k-kstart] = startx[2] + k*DKS(i,j)*dx[2];
    Gamma[k-kstart] = Gamma2s[lev][k];
    beta [k-kstart] = beta2s [lev][k];
    alpha[k-kstart] = alpha2s[lev][k];
  }

  return;
}


void init_volume()
{
  int ii,jj,kk;
  double x0lo,x0hi,x1lo,x1hi,r2dr_int,sthdth_int;

  //#pragma omp parallel for private(ii,jj,kk,x0lo,x0hi,r2dr_int,x1lo,x1hi,sthdth_int) schedule(guided)
  ISGLOOP(ii) {
    x0lo = startx[0] +  ii   *dx[0];
    x0hi = startx[0] + (ii+1)*dx[0];
    r2dr_int = romb(x0lo,x0hi,r2dr_dx);

    #if (NDIM>1)
    JSGLOOP(ii,jj) {
      x1lo = startx[1] +  jj   *DJS(ii)*dx[1];
      x1hi = startx[1] + (jj+1)*DJS(ii)*dx[1];
    #else
      x1lo = -1.0;
      x1hi =  1.0;
    #endif
      sthdth_int = romb(x1lo, x1hi, sthdth_dx);
      #if (NDIM==3)
      KSGLOOP(ii,jj,kk) {
      #endif
        ND_ELEM_LINEAR(geom,ii,jj,kk).volume = r2dr_int*sthdth_int*DKS(ii,jj)*dx[2];
      #if (NDIM==3)
      }
      #endif
    #if (NDIM>1)
    }
    #endif
  }
  return;
}

/* Initialize areas:  Area in the x^l-direction is defined as the integral of sqrt(g)
 *   over x^m and x^n.
 */
void init_area()
{
  int ii,jj,kk,kkk;
  double x0lo,x0hi,x1lo,x1hi,r2dr_int,sthdth_int,dph_int;

  //#pragma omp parallel for private(ii,jj,kk,x0lo,x0hi,x1lo,x1hi,r2dr_int,sthdth_int,dph_int) schedule(guided)
  for (ii=istart[0]; ii<=istop[0]; ii++) {
    x0lo = startx[0] +  ii   *dx[0];
    x0hi = startx[0] + (ii+1)*dx[0];
    r2dr_int = romb(x0lo,x0hi,r2dr_dx);

    #if (NDIM>1)
    for (jj=JS(ii,istart[1]); jj<=JS(ii,istop[1]); jj++) {
      x1lo = startx[1] +  jj   *DJS(ii)*dx[1];
      x1hi = startx[1] + (jj+1)*DJS(ii)*dx[1];
    #else
      x1lo = -1.0;
      x1hi =  1.0;
    #endif

      sthdth_int = romb(x1lo,x1hi,sthdth_dx);

      #if (NDIM==3)
      for (kk=KS(ii,jj,istart[2]); kk<=KS(ii,jj,istop[2]); kk++) {
      #endif

        dph_int = DKS(ii,jj)*dx[2];

        ND_ELEM_LINEAR(geom,ii,jj,kk).area[0] = r2dr_dx(x0lo)*sthdth_int*dph_int;

        #if (NDIM>1)
        ND_ELEM_LINEAR(geom,ii,jj,kk).area[1] = r2dr_int*sthdth_dx(x1lo)*dph_int;
        #endif

        #if (NDIM==3)
        ND_ELEM_LINEAR(geom,ii,jj,kk).area[2] = r2dr_int*sthdth_int;
        #endif

      #if (NDIM==3)
      }
      #endif
    #if (NDIM>1)
    }
    #endif
  }

  // make sure area on axis is a hard-zero
  #if (NDIM>1)
  if (startx[1] == -1 && istart[1] == 0) {
    for (ii=istart[0]; ii<=istop[0]; ii++) {
      jj = 0;
      #if (NDIM==3)
      for (kk=KS(ii,jj,istart[2]); kk<=KS(ii,jj,istop[2]); kk++) {
      #endif
        ND_ELEM_LINEAR(geom,ii,jj,kk).area[1] = 0.0;
      #if (NDIM==3)
      }
      #endif
    }
  }
  if (fabs(startx[1] + n2*dx[1] - 1) < 1.0e-15 && istop[1] == n2) {
    for (ii=istart[0]; ii<=istop[0]; ii++) {
      jj = JS(ii,istop[1]);
      #if (NDIM==3)
      for (kk=KS(ii,jj,istart[2]); kk<=KS(ii,jj,istop[2]); kk++) {
      #endif
        ND_ELEM_LINEAR(geom,ii,jj,kk).area[1] = 0.0;
      #if (NDIM==3)
      }
      #endif
    }
  }
  #endif

  return;
}

/* Initialize the (physical) volume-averaged connection coefficients (Christoffel symbols).
 *   Recall that there is no sum on repeated indices for Christoffel symbols.
 *   conn[l][m][n] = \Gamma^l_mn = { l }
 *                                 {m n}
 *   Christoffel symbols are symmetric in lower indices, i.e., conn[l][m][n] = conn[l][n][m].
 *   In this coordinate system, the metric elements are g_11(x1) = (dr/dx1)^2, g_22(x1,x2) = r^2 (dth/dx2)^2,
 *     and g_33(x1,x2) = r^2 sin^2(th), and g_ij=0 for i!=j, hence g_ij is orthogonal.
 *   Note that the indices of the C arrays are 0,1,2, but for the coordinates as written, they are 1,2,3.
 *   The quantities below are volume-averaged quantities, so first compute the Christoffel symbol, then
 *     take the volume average noting that:
 *     <f(x1,x2,x3)> = 1/dV \int f(x1,x2,x3) r^2(x1) sin(th(x2)) (dr/dx1)(x1) (dth/dx2)(x2) dx1 dx2 dx3
 *   For an orthogonal coordinate system, only the following symbols are non-zero (i != j):
 *     { i } = 0.5*(ln g_ii)_,i   ,   { i } = 0.5*(ln g_ii)_,j    ,   { i } = -0.5*(g_jj)_,i / g_ii
 *     {i i}                          {i j}                           {j j}
 */
void init_conn()
{
  int ii,jj,kk;
  double x0lo,x0hi,vrf,x1lo,x1hi,vthf;
  double rint_r2d2r_dx2,rint_r3,rint_rdr_dx2;
  double thint_sthdth_dx3,thint_sthd2th_dx2,thint_cthdth_dx2;
  int kp,nkp,k;
  double dA,dV;

  //#pragma omp parallel for private(ii,jj,kk,x0lo,x0hi,vrf,x1lo,x1hi,vthf,rint_r2d2r_dx2,rint_r3,rint_rdr_dx2,thint_sthdth_dx3,thint_sthd2th_dx2,thint_cthdth_dx2,kp,nkp,k,dA,dV) schedule(guided)
  ISLOOP(ii) {
    x0lo = startx[0] +  ii   *dx[0];
    x0hi = startx[0] + (ii+1)*dx[0];
    rint_r2d2r_dx2 = romb(x0lo,x0hi,r2d2r_dx2);
    rint_r3        = romb(x0lo,x0hi,r3       );
    rint_rdr_dx2   = romb(x0lo,x0hi,rdr_dx2  );
    vrf  = Vr(x0lo);
    #if (NDIM>1)
    JSLOOP(ii,jj) {
      x1lo = startx[1] +  jj   *DJS(ii)*dx[1];
      x1hi = startx[1] + (jj+1)*DJS(ii)*dx[1];
      #else
      x1lo = -1.0;
      x1hi =  1.0;
      #endif
      vthf              = romb(x1lo,x1hi,sthdth_dx  );
      thint_sthdth_dx3  = romb(x1lo,x1hi,sthdth_dx3 );
      thint_sthd2th_dx2 = romb(x1lo,x1hi,sthd2th_dx2);
      thint_cthdth_dx2  = romb(x1lo,x1hi,cthdth_dx2 );
      #if (NDIM==3)
      KSLOOP(ii,jj,kk) {
      #endif

        memset(&ND_ELEM_LINEAR(geom,ii,jj,kk).conn[0][0][0], 0, 27*sizeof(double));

        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[0][0][0] = rint_r2d2r_dx2/vrf;
        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[0][1][1] = -(rint_r3/vrf)*(thint_sthdth_dx3/vthf);

        // This angular integral is analytic for all th(x), but who cares?
        //ND_ELEM_LINEAR(geom,ii,jj,kk).conn[0][2][2] = -(int_r3(xr0)/vrf)*(romb(xt0,xt0+spider_fact[ii]*dx[1],sth3dth_dx)/vthf);

        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[1][0][1] = rint_rdr_dx2/vrf;
        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[1][1][0] = ND_ELEM_LINEAR(geom,ii,jj,kk).conn[1][0][1];
        // A.S.:  Should we force [1][1][0] or [0][1][1] to give 0 as in the [0][2][2] and [1][2][2] cases below?

        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[2][0][2] = ND_ELEM_LINEAR(geom,ii,jj,kk).conn[1][0][1];
        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[2][2][0] = ND_ELEM_LINEAR(geom,ii,jj,kk).conn[2][0][2];

        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[2][1][2] = thint_cthdth_dx2/vthf;
        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[2][2][1] = ND_ELEM_LINEAR(geom,ii,jj,kk).conn[2][1][2];

        // Set certain coefficients to ensure exact angular momentum conservation.
        // For the phi-momentum equation (m=3), we want the geometric source term to vanish in order
        // to conserve phi-angular momentum.  That is, we want
        // 0 = Gamma^l_n3 T^n_l
        //   = Gamma^3_13*rho*v^1*v_3 + Gamma^3_23*rho*v^2*v_3 + Gamma^1_33*rho*v^3*v_1 + Gamma^2_33*rho*v^3*v_2,
        //   = Gamma^3_13*v^1*g_33*v^3 + Gamma^3_23*v^2*g_33*v^3 + Gamma^1_33*v^3*g_11*v^1 + Gamma^2_33*v^3*g_22*v^2,
        //   = (Gamma^3_13*g_33 + Gamma^1_33*g_11)*v^1 + (Gamma^3_23*g_33 + Gamma^2_33*g_22)*v^2
        // Since v^1 and v^2 are arbitrary, it must be that the expressions in ()'s are each = 0.
        // Thus, we can set Gamma^1_33 = -(g_33/g_11)*Gamma^3_13 and Gamma^2_33 = -(g_33/g_22)*Gamma^3_23.
        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[0][2][2] = -ND_ELEM_LINEAR(geom,ii,jj,kk).conn[2][2][0] * ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][2]/ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][0];
        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[1][2][2] = -ND_ELEM_LINEAR(geom,ii,jj,kk).conn[2][2][1] * ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][2]/ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][1];

        // Set Gamma^2_22 to get (\grad P)_\theta = 0 correct when it should be
        #if (NDIM>1)
        if (DKS(ii,jj+1) < DKS(ii,jj)) {
          // Refinement boundary, 2 outer neighbors, shifted index, original area
          dA = ND_ELEM_LINEAR(geom,ii,jj+1,2*kk).area[1] + ND_ELEM_LINEAR(geom,ii,jj+1,2*kk+1).area[1];
        } else if (DKS(ii,jj+1) > DKS(ii,jj)) {
          // Coarsening boundary, 1 outer neighbor, shifted index, half area
          dA = 0.5*ND_ELEM_LINEAR(geom,ii,jj+1,kk/2).area[1];
        } else {
          // Regular boundary, 1 outer neighbor, original index, original area
          dA = ND_ELEM_LINEAR(geom,ii,jj+1,kk).area[1];
        }
        dA -= ND_ELEM_LINEAR(geom,ii,jj,kk).area[1];
        dV  = ND_ELEM_LINEAR(geom,ii,jj,kk).volume;
        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[1][1][1] = dA/dV - ND_ELEM_LINEAR(geom,ii,jj,kk).conn[2][1][2];

        #else
        ND_ELEM_LINEAR(geom,ii,jj,kk).conn[1][1][1] = thint_sthd2th_dx2/vthf;
        #endif

      #if (NDIM==3)
      }
      #endif
    #if (NDIM>1)
    }
    #endif
  }

  return;
}


/* Initialize the covariant components, g_ij, of the metric tensor.  Since we implicitly
 * assume our coordinate system is orthogonal, it follows that g_ij is a diagonal tensor,
 * i.e., g_ij=0 when i!=j.  Thus, we only need to store the diagonal components, g_ii,
 * for i=0,1,2.  We need zone volume averages (gcov[0][i]) and interface area averages
 * (gcov[k+1][i] for k=0,1,2 representing the direction of the interface normal) of g_ii.
 * The components are:  g_00=[dr/dx0]^2, g_11=[r*dth/dx1]^2, g_22=[r*sin(th)]^2.
 * For example, to get gcov[0][0], which is the zone volume average of g_00, take
 * <g_00> = \int_x0lo^x0hi [dr/dx0]^2 r^2 dr/dx0 dx0 / \int_x0lo^x0hi r^2 dr/dx0 dx0.
 */
void init_gcov()
{
  double x0lo,x0hi,vrf,x1lo,x1hi,vthf;
  double rint_r2dr_dx3,thint_sthdth_dx3,thint_sth3dth_dx;
  int ii,jj,kk,loc,m;

  //#pragma omp parallel for private(ii,jj,kk,loc,m,x0lo,x0hi,vrf,x1lo,x1hi,vthf,rint_r2dr_dx3,thint_sthdth_dx3,thint_sth3dth_dx) schedule(guided)
  ISGLOOP(ii) {
    x0lo = startx[0] +  ii   *dx[0];
    x0hi = startx[0] + (ii+1)*dx[0];
    vrf           = romb(x0lo,x0hi,r2dr_dx );
    rint_r2dr_dx3 = romb(x0lo,x0hi,r2dr_dx3);

    #if (NDIM>1)
    JSGLOOP(ii,jj) {
      x1lo = startx[1] +  jj   *DJS(ii)*dx[1];
      x1hi = startx[1] + (jj+1)*DJS(ii)*dx[1];
    #else
      x1lo = -1.0;
      x1hi =  1.0;
    #endif
      vthf             = romb(x1lo,x1hi,sthdth_dx );
      thint_sthdth_dx3 = romb(x1lo,x1hi,sthdth_dx3);
      thint_sth3dth_dx = romb(x1lo,x1hi,sth3dth_dx);

      #if (NDIM==3)
      KSGLOOP(ii,jj,kk) {
      #endif

        // initialize to zero
        for (loc=0; loc<=NDIM; loc++) {
          for (m=0; m<SPACEDIM; m++) {
            ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[loc][m] = 0.0;
          }
        }

        // volume-averaged
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][0] = rint_r2dr_dx3/vrf;
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][1] = fabs(0.2*(pow(r_of_x(rx_info,x0hi),5) - pow(r_of_x(rx_info,x0lo),5))/vrf * thint_sthdth_dx3/vthf);
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][2] = fabs(0.2*(pow(r_of_x(rx_info,x0hi),5) - pow(r_of_x(rx_info,x0lo),5))/vrf * thint_sth3dth_dx/vthf);

        // 0-direction -- face-averaged
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[1][0] = pow(dr_dx(x0lo),2);
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[1][1] = pow(r_of_x(rx_info,x0lo),2)*thint_sthdth_dx3/vthf;
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[1][2] = pow(r_of_x(rx_info,x0lo),2)*thint_sth3dth_dx/vthf;

        #if (NDIM>1)
        // 1-direction -- face-averaged
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[2][0] = ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][0];
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[2][1] = fabs(0.2*(pow(r_of_x(rx_info,x0hi),5) - pow(r_of_x(rx_info,x0lo),5))/vrf * pow(dth_dx(x1lo),2));

        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[2][2] = fabs(0.2*(pow(r_of_x(rx_info,x0hi),5) - pow(r_of_x(rx_info,x0lo),5))/vrf * pow(sin(th_of_x(thx_info,x1lo)),2));
        #endif

        #if (NDIM==3)
        // 2-direction -- face-averaged
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[3][0] = ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][0];
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[3][1] = ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][1];
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[3][2] = ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][2];
        #endif
      #if (NDIM==3)
      }
      #endif
    #if (NDIM>1)
    }
    #endif
  }

  return;
}

/* Initialize the contravariant components of g.  For a diagonal metric, these are
 * just the inverses of the covariant components, previously computed above.  We only
 * need zone volume averages (loc=0) of these components.
 */
void init_gcon()
{
  int ii,jj,kk,i;

  //#pragma omp parallel for private(ii,jj,kk,i)
  ZGLOOP {
    for (i=0; i<SPACEDIM; i++) {
      if (fabs(ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][i]) > 1.0e-50) {
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcon[i] = 1.0/ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][i];
      } else {
        ND_ELEM_LINEAR(geom,ii,jj,kk).gcon[i] = 0.0;
      }
    }
  }

  return;
}

/* Initialize the metric scale factors, defined as sqrt(g_ii), for i=0,1,2.
 * These are needed as zone volume averages (loc=0) and interface area averages (loc=1,2,3).
 */
void init_scale()
{
  int ii,jj,kk,loc,dd;

  //#pragma omp parallel for private(ii,jj,kk,loc,dd)
  ZGLOOP {
    for (loc=0; loc<=NDIM; loc++) {
      SLOOP {
        ND_ELEM_LINEAR(geom,ii,jj,kk).scale[loc][dd] = sqrt(ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[loc][dd]);
      }
    }
  }

  return;
}
