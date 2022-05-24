
#include "decs.h"

inline void p_uflux(const int nhydro,const double* restrict p, const zone_geom* restrict g, double press, const int dir, double* restrict u, double* restrict f)
{
	int dd,vv;
	double vcov[SPACEDIM],vsq,etot;
	double v = p[U1+dir];

	geom_lower(&p[U1], vcov, g->gcov[1+dir]);

	u[RHO] = p[RHO];
	vsq = 0.0;
	for (dd=0; dd<SPACEDIM; dd++) {
		vsq += p[U1+dd]*vcov[dd];
	}
	etot = 0.5*u[RHO]*vsq + p[UU];
	u[ETOT] = etot;
	for (dd=U1; dd<U1+SPACEDIM; dd++) {
		u[dd] = u[RHO]*vcov[dd-U1];
	}
	for (vv=U1+SPACEDIM; vv<nhydro; vv++) {
		u[vv] = u[RHO]*p[vv];
	}

	f[RHO] = u[RHO]*v;
	f[ETOT] = v*(etot + press);
	for (dd=U1; dd<U1+SPACEDIM; dd++) {
		f[dd] = f[RHO]*vcov[dd-U1];
	}
	f[U1+dir] += press;
	for (vv=U1+SPACEDIM; vv<nhydro; vv++) {
		f[vv] = f[RHO]*p[vv];
	}

	return;
}

inline void pflux(const int nhydro,const double* restrict p, double press, const zone_geom* restrict g, const int dir, double* restrict f)
{
	int dd,vv;
	double vcov[SPACEDIM],vsq,etot;
	double v = p[U1+dir];
	double area = g->area[dir];

	geom_lower(&p[U1], vcov, g->gcov[1+dir]);

	f[RHO] = area*p[RHO]*v;
	for (dd=U1; dd<U1+SPACEDIM; dd++) {
		f[dd] = f[RHO]*vcov[dd-U1];
	}
	f[U1+dir] += area*press;
	vsq = 0.0;
	for (dd=0; dd<SPACEDIM; dd++) vsq += p[U1+dd]*vcov[dd];
	etot = 0.5*p[RHO]*vsq + p[UU];
	f[ETOT] = area*v*(etot + press);
	for (vv=U1+SPACEDIM; vv<nhydro; vv++) {
		f[vv] = f[RHO]*p[vv];
	}

	return;
}

//int isGlobal(void *p) { return 0; }
//#pragma omp begin declare variant match(device={kind(gpu)})
//int isGlobal(void *p) { return __nvvm_isspace_global(p); }
//#pragma omp end declare variant

inline int riemann_flux_LLF(int nhydro, int nvars,double pl[], double pr[], const zone_geom const *g, const int dir, double flux[], double *sig_speed, double *vriemann)
{
	int vv,dd;
	double ul[nhydro],ur[nhydro],fl[nhydro],fr[nhydro];
	double vsql, vsqr;
	double p_l, p_r, cs_l, cs_r;
	double Gamma_l, Gamma_r;
	double u_l,u_r;
	double cmaxl,cmaxr,cmax;
	double area = g->area[dir];

	if (area < 1.0e-16) {
		*sig_speed = 1.0e-200;
		SLOOP {
			vriemann[dd] = 0.0;
		}
		HLOOP {
			flux[vv] = 0.0;
		}
		return 0;
	}

    //TIMER_START("riemann_flux_LLF");

	vv = nvars;
	p_l = pl[vv];		p_r = pr[vv];	vv++;
	Gamma_l = pl[vv];	Gamma_r = pr[vv];

	cs_l = sqrt(Gamma_l*p_l/pl[RHO]);
	cs_r = sqrt(Gamma_r*p_r/pr[RHO]);

	cs_l /= g->scale[1+dir][dir];
	cs_r /= g->scale[1+dir][dir];

        //nhydro p g press dir u f
	//p_uflux(nhydro,pl, g, p_l, dir, ul, fl);
	//p_uflux(nhydro,pr, g, p_r, dir, ur, fr);

	  double vcov[SPACEDIM],vsq,etot;
	  double v = pl[U1+dir];
	  geom_lower(&pl[U1], vcov, g->gcov[1+dir]);
	  ul[RHO] = pl[RHO];
	  vsq = 0.0;
	  for (dd=0; dd<SPACEDIM; dd++) {
	  	vsq += pl[U1+dd]*vcov[dd];
	  }
	  etot = 0.5*ul[RHO]*vsq + pl[UU];
	  ul[ETOT] = etot;
	  for (dd=U1; dd<U1+SPACEDIM; dd++) {
	  	ul[dd] = ul[RHO]*vcov[dd-U1];
	  }
	  for (vv=U1+SPACEDIM; vv<nhydro; vv++) {
	  	ul[vv] = ul[RHO]*pl[vv];
	  }
	  fl[RHO] = ul[RHO]*v;
	  fl[ETOT] = v*(etot + p_l);
	  for (dd=U1; dd<U1+SPACEDIM; dd++) {
	  	fl[dd] = fl[RHO]*vcov[dd-U1];
	  }
	  fl[U1+dir] += p_l;
	  for (vv=U1+SPACEDIM; vv<nhydro; vv++) {
	  	fl[vv] = fl[RHO]*pl[vv];
	  }

	  v = pl[U1+dir];
	  geom_lower(&pr[U1], vcov, g->gcov[1+dir]);
	  ur[RHO] = pr[RHO];
	  vsq = 0.0;
	  for (dd=0; dd<SPACEDIM; dd++) {
	  	vsq += pr[U1+dd]*vcov[dd];
	  }
	  etot = 0.5*ur[RHO]*vsq + pr[UU];
	  ur[ETOT] = etot;
	  for (dd=U1; dd<U1+SPACEDIM; dd++) {
	  	ur[dd] = ur[RHO]*vcov[dd-U1];
	  }
	  for (vv=U1+SPACEDIM; vv<nhydro; vv++) {
	  	ur[vv] = ur[RHO]*pr[vv];
	  }
	  fr[RHO] = ur[RHO]*v;
	  fr[ETOT] = v*(etot + p_r);
	  for (dd=U1; dd<U1+SPACEDIM; dd++) {
	  	fr[dd] = fr[RHO]*vcov[dd-U1];
	  }
	  fr[U1+dir] += p_r;
	  for (vv=U1+SPACEDIM; vv<nhydro; vv++) {
	  	fr[vv] = fr[RHO]*pr[vv];
	  }

	u_l = pl[U1+dir];
	u_r = pr[U1+dir];

	cmaxl = MAX(fabs(u_l+cs_l), fabs(u_l-cs_l));
	cmaxr = MAX(fabs(u_r+cs_r), fabs(u_r-cs_r));
	cmax = MAX(cmaxl, cmaxr);
	*sig_speed = cmax;

	SLOOP {
		vriemann[dd] = 0.5*(pl[U1+dd] + pr[U1+dd]);
	}

	for(vv=0;vv<nhydro;vv++) {
                flux[vv] = 0.5*area*(fl[vv] + fr[vv] - cmax*(ur[vv] - ul[vv]));
	}

    //TIMER_STOP;
	return 0;
}

int riemann_flux_HLLE(int nhydro,int nvars,double pl[], double pr[], const zone_geom *g, const int dir, double flux[], double *sig_speed, double *vriemann)
{
  int vv,dd;
  double ul[NSIZE],ur[NSIZE],fl[NSIZE],fr[NSIZE],du[NSIZE];
  double vsql,vsqr;
  double p_l,p_r,cs_l,cs_r,p_s,plr;
  double p_max,p_min,pratio;
  double Gamma_l,Gamma_r,Gamma_avg;
  double a_l,a_r,b_l,b_r,g_l,g_r;
  double u_l,u_r,q_l,q_r,sl,sr,ss;
  double cmaxl,cmaxr,cmax,cmin;
  double area = g->area[dir];
  double rhowl,rhowr;

  //TIMER_START("riemann_flux_HLLE");

  vv = nvars;
  p_l = pl[vv];  p_r = pr[vv];  vv++;
  Gamma_l = pl[vv];  Gamma_r = pr[vv];

  if (area < 1.0e-16) {
    *sig_speed = 1.0e-200;
    SLOOP {
      vriemann[dd] = 0.5*(pl[U1+dd] + pr[U1+dd]);
    }
    HLOOP {
      flux[vv] = 0.0;
    }

    return 0;
  }

  u_l = pl[U1+dir]*g->scale[1+dir][dir];
  u_r = pr[U1+dir]*g->scale[1+dir][dir];

  cs_l = sqrt(Gamma_l*p_l/pl[RHO]);
  cs_r = sqrt(Gamma_r*p_r/pr[RHO]);

  if (p_l > p_r) {
    p_max = p_l;
    p_min = p_r;
  } else {
    p_max = p_r;
    p_min = p_l;
  }
  pratio = p_max/p_min;

  p_s = 0.5*(p_l + p_r) - 0.125*(u_r-u_l)*(pl[RHO]+pr[RHO])*(cs_l+cs_r);
  p_s = MAX(p_s, 0.0);

  // if (pratio > 2.0) {
  //   if (p_s < p_min) {
  //     // two rarefaction waves
  //     double Gamma_avg = 0.5*(Gamma_l + Gamma_r);
  //     double plr = pow(p_l/p_r, (Gamma_avg-1.0)/(2.0*Gamma_avg));
  //     ss = (plr*u_l/cs_l + u_r/cs_r + 2.0*(plr - 1.0)/(Gamma_avg - 1.0))/(plr/cs_l + 1.0/cs_r);
  //     p_s = 0.5*(p_l*pow(1.0 + (Gamma_l-1)/(2.0*cs_l)*(u_l - ss),2*Gamma_l/(Gamma_l-1.0)) + p_r*pow(1.0 + (Gamma_r-1.0)/(2.0*cs_r)*(ss - u_r),2.0*Gamma_r/(Gamma_r-1.0)));
  //     p_s = MAX(0.0, p_s);
  //     q_l = (p_s <= p_l) ? 1 : sqrt(1 + (Gamma_l+1)/(2*Gamma_l)*(p_s/p_l - 1.0));
  //     q_r = (p_s <= p_r) ? 1 : sqrt(1 + (Gamma_r+1)/(2*Gamma_r)*(p_s/p_r - 1.0));
  //     sl = u_l - cs_l*q_l;
  //     sr = u_r + cs_r*q_r;
  //   } else {
  //     // two shock waves
  //     double a_l = 2.0/((Gamma_l+1.0)*pl[RHO]);
  //     double a_r = 2.0/((Gamma_r+1.0)*pr[RHO]);
  //     double b_l = (Gamma_l-1.0)/(Gamma_l+1.0)*p_l;
  //     double b_r = (Gamma_r-1.0)/(Gamma_r+1.0)*p_r;
  //     double g_l = sqrt(a_l/(p_s + b_l));
  //     double g_r = sqrt(a_r/(p_s + b_r));
  //     p_s = (g_l*p_l + g_r*p_r - (u_r-u_l))/(g_l + g_r);
  //     p_s = MAX(0.0, p_s);
  //     q_l = (p_s <= p_l) ? 1 : sqrt(1 + (Gamma_l+1)/(2*Gamma_l)*(p_s/p_l - 1.0));
  //     q_r = (p_s <= p_r) ? 1 : sqrt(1 + (Gamma_r+1)/(2*Gamma_r)*(p_s/p_r - 1.0));
  //     sl = u_l - cs_l*q_l;
  //     sr = u_r + cs_r*q_r;
  //   }
  // } else {
  //   q_l = (p_s <= p_l) ? 1 : sqrt(1 + (Gamma_l+1)/(2*Gamma_l)*(p_s/p_l - 1.0));
  //   q_r = (p_s <= p_r) ? 1 : sqrt(1 + (Gamma_r+1)/(2*Gamma_r)*(p_s/p_r - 1.0));
  //   sl = u_l - cs_l*q_l;
  //   sr = u_r + cs_r*q_r;
  // }

  if (pratio < 2.0 && p_s >= p_min && p_s <= p_max) {
    q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
    q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
    // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
    // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
    sl = u_l - cs_l*q_l;
    sr = u_r + cs_r*q_r;
    ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
  } else {
    if (p_s < p_min) {
      // Two rarefaction waves (see Toro, eqs. 9.35-9.36)
      Gamma_avg = 0.5*(Gamma_l + Gamma_r);
      plr = pow(p_l/p_r, (Gamma_avg-1.0)/(2.0*Gamma_avg));
      ss = (plr*u_l/cs_l + u_r/cs_r + 2.0*(plr - 1.0)/(Gamma_avg - 1.0))/(plr/cs_l + 1.0/cs_r);
      p_s = 0.5*(p_l*pow(1.0 + (Gamma_l-1.0)/(2.0*cs_l)*(u_l - ss),2.0*Gamma_l/(Gamma_l-1.0)) + p_r*pow(1.0 + (Gamma_r-1.0)/(2.0*cs_r)*(ss - u_r),2.0*Gamma_r/(Gamma_r-1.0)));
      p_s = MAX(p_s, 0.0);
      q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
      q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
      // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
      // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
      sl = u_l - cs_l*q_l;
      sr = u_r + cs_r*q_r;
      ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
    } else {
      // Two shock waves (see Toro, eqs. 9.42-9.43)
      a_l = 2.0/((Gamma_l+1.0)*pl[RHO]);
      a_r = 2.0/((Gamma_r+1.0)*pr[RHO]);
      b_l = (Gamma_l-1.0)/(Gamma_l+1.0)*p_l;
      b_r = (Gamma_r-1.0)/(Gamma_r+1.0)*p_r;
      g_l = sqrt(a_l/(p_s + b_l));
      g_r = sqrt(a_r/(p_s + b_r));
      p_s = (g_l*p_l + g_r*p_r - (u_r-u_l))/(g_l + g_r);
      p_s = MAX(p_s, 0.0);
      q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
      q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
      // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
      // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
      sl = u_l - cs_l*q_l;
      sr = u_r + cs_r*q_r;
      ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
    }
  }
  
  rhowl = sqrt(pl[RHO]);
  rhowr = sqrt(pr[RHO]);
  SLOOP {
    vriemann[dd] = (rhowl*pl[U1+dd] + rhowr*pr[U1+dd])/(rhowl+rhowr);
  }

  p_uflux(nhydro,pl,g,p_l,dir,ul,fl);
  p_uflux(nhydro,pr,g,p_r,dir,ur,fr);

  cmax = MAX(sr/g->scale[1+dir][dir], 0.0);
  cmin = fabs(MIN(sl/g->scale[1+dir][dir], 0.0));
  *sig_speed = MAX(cmax,cmin);

  HLOOP {
    flux[vv] = area*(cmin*fr[vv] + cmax*fl[vv] - cmax*cmin*(ur[vv]-ul[vv]))/(cmax+cmin);
  }

  //TIMER_STOP;

  return 0;
}


inline int riemann_flux_HLLC(const int nhydro, const int nvars,const double* restrict pl, const double* restrict pr, const zone_geom* restrict g, const int dir, 
                      double* restrict flux, double* restrict sig_speed, double* restrict vriemann)
{
  int vv,dd,ssflag;
  double u_l,u_r;
  double p_l,p_r,cs_l,cs_r,cmax,plr;
  double Gamma_l,Gamma_r;
  double temp_l,temp_r;
  double ent_l,ent_r;
  double sl,sr,ss,p_s,slr;
  double q_l,q_r,pratio,p_max,p_min;
  double a_l,a_r,b_l,b_r,g_l,g_r,z;
  double vsq,Gamma_max,Gamma_min,Gamma_avg;
  double u[NSIZE];
  double fx[NSIZE];
  double ds[NSIZE];
  double Ul[NSIZE],Ur[NSIZE],fl[NSIZE],fr[NSIZE];
  double Us[NSIZE];
  double vcov[3];
  double area = g->area[dir];
  double eos[NEOS];
  double uint,pnew;
  double ustar;
  double sonicl,sonicr;
  int isshock = 0;

  //TIMER_START("riemann_flux_HLLC");

  vv = nvars;
  p_l = pl[vv];  p_r = pr[vv];  vv++;
  Gamma_l = pl[vv];  Gamma_r = pr[vv];

  if (area < 1.0e-16) {
    *sig_speed = 1.0e-200;
    SLOOP { vriemann[dd] = 0.5*(pl[U1+dd] + pr[U1+dd]); }
    HLOOP { flux[vv] = 0.0; }
    //TIMER_STOP;
    return isshock; 
  }

  u_l = pl[U1+dir]*g->scale[1+dir][dir];
  u_r = pr[U1+dir]*g->scale[1+dir][dir];

  cs_l = sqrt(Gamma_l*p_l/pl[RHO]);
  cs_r = sqrt(Gamma_r*p_r/pr[RHO]);

  if (p_l > p_r) {
    p_max = p_l;
    p_min = p_r;
  } else {
    p_max = p_r;
    p_min = p_l;
  }
  pratio = p_max/p_min;

  isshock = ((p_max-p_min)/p_min > 0.5) ? 1 : 0;

  p_s = 0.5*(p_l + p_r) - 0.125*(u_r-u_l)*(pl[RHO]+pr[RHO])*(cs_l+cs_r);
  p_s = MAX(0,p_s);

  if (pratio < 2.0 && p_s >= p_min && p_s <= p_max) {
    q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
    q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
    // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
    // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
    sl = u_l - cs_l*q_l;
    sr = u_r + cs_r*q_r;
    ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
  } else {
    if (p_s < p_min) {
      // Two rarefaction waves (see Toro, eqs. 9.35-9.36)
      Gamma_avg = 0.5*(Gamma_l + Gamma_r);
      plr = pow(p_l/p_r, (Gamma_avg-1.0)/(2.0*Gamma_avg));
      ss = (plr*u_l/cs_l + u_r/cs_r + 2.0*(plr - 1.0)/(Gamma_avg - 1.0))/(plr/cs_l + 1.0/cs_r);
      p_s = 0.5*(p_l*pow(1.0 + (Gamma_l-1.0)/(2.0*cs_l)*(u_l - ss),2.0*Gamma_l/(Gamma_l-1.0)) + p_r*pow(1.0 + (Gamma_r-1.0)/(2.0*cs_r)*(ss - u_r),2.0*Gamma_r/(Gamma_r-1.0)));
      p_s = MAX(p_s, 0.0);
      q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
      q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
      // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
      // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
      sl = u_l - cs_l*q_l;
      sr = u_r + cs_r*q_r;
      ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
    } else {
      // Two shock waves (see Toro, eqs. 9.42-9.43)
      a_l = 2.0/((Gamma_l+1.0)*pl[RHO]);
      a_r = 2.0/((Gamma_r+1.0)*pr[RHO]);
      b_l = (Gamma_l-1.0)/(Gamma_l+1.0)*p_l;
      b_r = (Gamma_r-1.0)/(Gamma_r+1.0)*p_r;
      g_l = sqrt(a_l/(p_s + b_l));
      g_r = sqrt(a_r/(p_s + b_r));
      p_s = (g_l*p_l + g_r*p_r - (u_r-u_l))/(g_l + g_r);
      p_s = MAX(p_s, 0.0);
      q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
      q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
      // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
      // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
      sl = u_l - cs_l*q_l;
      sr = u_r + cs_r*q_r;
      ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
    }
  }

  if (ss > sr) sr = 1.001*ss;
  if (ss < sl) sl = 1.001*ss;

  *sig_speed = MAX(fabs(sl),fabs(sr))/g->scale[1+dir][dir];

  // supersonic cases
  if (sl > 0.0) {
    pflux(nhydro,pl, p_l, g, dir, flux);
    SLOOP {
      vriemann[dd] = pl[U1+dd];
    }
    //TIMER_STOP;
    return isshock;
  } else if (sr < 0.0) {
    pflux(nhydro,pr, p_r, g, dir, flux);
    SLOOP {
      vriemann[dd] = pr[U1+dd];
    }
    //TIMER_STOP;
    return isshock;
  }

  if (ss > 0.0) {
    p_uflux(nhydro,pl, g, p_l, dir, u, fx);
    Us[0] = pl[RHO]*(sl - u_l)/(sl - ss);
    Us[ETOT] = (u[ETOT]*(sl - u_l) - p_l*u_l + p_s*ss)/(sl - ss);
    SLOOP {
      Us[U1+dd] = pl[U1+dd];
    }
    Us[U1+dir] = ss/g->scale[1+dir][dir];
    geom_lower(&Us[U1], vcov, g->gcov[1+dir]);
    vsq = 0.0;
    SLOOP {
      vsq += Us[U1+dd]*vcov[dd];
      vriemann[dd] = Us[U1+dd];
      Us[U1+dd] = vcov[dd]*Us[0];
    }
    for (vv=U1+SPACEDIM;vv<nhydro;vv++) {
      Us[vv] = Us[0]*pl[vv];
    }
    slr = sl/g->scale[1+dir][dir];
  } else {
    p_uflux(nhydro,pr, g, p_r, dir, u, fx);
    Us[0] = pr[RHO]*(sr - u_r)/(sr - ss);
    Us[ETOT] = (u[ETOT]*(u_r - sr) + p_r*u_r - p_s*ss)/(ss - sr);
    SLOOP {
      Us[U1+dd] = pr[U1+dd];
    }
    Us[U1+dir] = ss/g->scale[1+dir][dir];
    geom_lower(&Us[U1], vcov, g->gcov[1+dir]);
    vsq = 0.0;
    SLOOP {
      vsq += Us[U1+dd]*vcov[dd];
      vriemann[dd] = Us[U1+dd];
      Us[U1+dd] = vcov[dd]*Us[0];
    }
    for (vv=U1+SPACEDIM;vv<nhydro;vv++) {
      Us[vv] = Us[0]*pr[vv];
    }
    slr = sr/g->scale[1+dir][dir];
  }

  HLOOP {
    flux[vv] = area*(fx[vv] + slr*(Us[vv] - u[vv]));
  }
  
  //TIMER_STOP;
  return isshock;
}



/*
void p_uflux(const int nhydro,const double* restrict p, const zone_geom* restrict g, double press, const int dir, double* restrict u, double* restrict f)
{
	int dd,vv;
	double vcov[SPACEDIM],vsq,etot;
	double v = p[U1+dir];

	geom_lower(&p[U1], vcov, g->gcov[1+dir]);

	u[RHO] = p[RHO];
	vsq = 0.0;
	for (dd=0; dd<SPACEDIM; dd++) {
		vsq += p[U1+dd]*vcov[dd];
	}
	etot = 0.5*u[RHO]*vsq + p[UU];
	u[ETOT] = etot;
	for (dd=U1; dd<U1+SPACEDIM; dd++) {
		u[dd] = u[RHO]*vcov[dd-U1];
	}
	for (vv=U1+SPACEDIM; vv<nhydro; vv++) {
		u[vv] = u[RHO]*p[vv];
	}

	f[RHO] = u[RHO]*v;
	f[ETOT] = v*(etot + press);
	for (dd=U1; dd<U1+SPACEDIM; dd++) {
		f[dd] = f[RHO]*vcov[dd-U1];
	}
	f[U1+dir] += press;
	for (vv=U1+SPACEDIM; vv<nhydro; vv++) {
		f[vv] = f[RHO]*p[vv];
	}

	return;
}

void pflux(const int nhydro,const double* restrict p, double press, const zone_geom* restrict g, const int dir, double* restrict f)
{
	int dd,vv;
	double vcov[SPACEDIM],vsq,etot;
	double v = p[U1+dir];
	double area = g->area[dir];

	geom_lower(&p[U1], vcov, g->gcov[1+dir]);

	f[RHO] = area*p[RHO]*v;
	for (dd=U1; dd<U1+SPACEDIM; dd++) {
		f[dd] = f[RHO]*vcov[dd-U1];
	}
	f[U1+dir] += area*press;
	vsq = 0.0;
	for (dd=0; dd<SPACEDIM; dd++) vsq += p[U1+dd]*vcov[dd];
	etot = 0.5*p[RHO]*vsq + p[UU];
	f[ETOT] = area*v*(etot + press);
	for (vv=U1+SPACEDIM; vv<nhydro; vv++) {
		f[vv] = f[RHO]*p[vv];
	}

	return;
}

//int isGlobal(void *p) { return 0; }
//#pragma omp begin declare variant match(device={kind(gpu)})
//int isGlobal(void *p) { return __nvvm_isspace_global(p); }
//#pragma omp end declare variant


int riemann_flux_LLF(int nhydro, int nvars,double pl[], double pr[], const zone_geom *g, const int dir, double flux[], double *sig_speed, double *vriemann)
{
	int vv,dd;
	double ul[nhydro],ur[nhydro],fl[nhydro],fr[nhydro];
	double vsql, vsqr;
	double p_l, p_r, cs_l, cs_r;
	double Gamma_l, Gamma_r;
	double u_l,u_r;
	double cmaxl,cmaxr,cmax;
	double area = g->area[dir];

	if (area < 1.0e-16) {
		*sig_speed = 1.0e-200;
		SLOOP {
			vriemann[dd] = 0.0;
		}
		HLOOP {
			flux[vv] = 0.0;
		}
		return 0;
	}

    //TIMER_START("riemann_flux_LLF");

	vv = nvars;
	p_l = pl[vv];		p_r = pr[vv];	vv++;
	Gamma_l = pl[vv];	Gamma_r = pr[vv];

	cs_l = sqrt(Gamma_l*p_l/pl[RHO]);
	cs_r = sqrt(Gamma_r*p_r/pr[RHO]);

	cs_l /= g->scale[1+dir][dir];
	cs_r /= g->scale[1+dir][dir];

	p_uflux(nhydro,pl, g, p_l, dir, ul, fl);
	p_uflux(nhydro,pr, g, p_r, dir, ur, fr);

	u_l = pl[U1+dir];
	u_r = pr[U1+dir];

	cmaxl = MAX(fabs(u_l+cs_l), fabs(u_l-cs_l));
	cmaxr = MAX(fabs(u_r+cs_r), fabs(u_r-cs_r));
	cmax = MAX(cmaxl, cmaxr);
	*sig_speed = cmax;

	SLOOP {
		vriemann[dd] = 0.5*(pl[U1+dd] + pr[U1+dd]);
	}

	for(vv=0;vv<nhydro;vv++) {
                flux[vv] = 0.5*area*(fl[vv] + fr[vv] - cmax*(ur[vv] - ul[vv]));
	}

    //TIMER_STOP;
	return 0;
}

int riemann_flux_HLLE(double pl[], double pr[], const zone_geom *g, const int dir, double flux[], double *sig_speed, double *vriemann)
{
  int vv,dd;
  double ul[nhydro],ur[nhydro],fl[nhydro],fr[nhydro],du[nhydro];
  double vsql,vsqr;
  double p_l,p_r,cs_l,cs_r,p_s,plr;
  double p_max,p_min,pratio;
  double Gamma_l,Gamma_r,Gamma_avg;
  double a_l,a_r,b_l,b_r,g_l,g_r;
  double u_l,u_r,q_l,q_r,sl,sr,ss;
  double cmaxl,cmaxr,cmax,cmin;
  double area = g->area[dir];
  double rhowl,rhowr;

  TIMER_START("riemann_flux_HLLE");

  vv = nvars;
  p_l = pl[vv];  p_r = pr[vv];  vv++;
  Gamma_l = pl[vv];  Gamma_r = pr[vv];

  if (area < 1.0e-16) {
    *sig_speed = 1.0e-200;
    SLOOP {
      vriemann[dd] = 0.5*(pl[U1+dd] + pr[U1+dd]);
    }
    HLOOP {
      flux[vv] = 0.0;
    }

    return 0;
  }

  u_l = pl[U1+dir]*g->scale[1+dir][dir];
  u_r = pr[U1+dir]*g->scale[1+dir][dir];

  cs_l = sqrt(Gamma_l*p_l/pl[RHO]);
  cs_r = sqrt(Gamma_r*p_r/pr[RHO]);

  if (p_l > p_r) {
    p_max = p_l;
    p_min = p_r;
  } else {
    p_max = p_r;
    p_min = p_l;
  }
  pratio = p_max/p_min;

  p_s = 0.5*(p_l + p_r) - 0.125*(u_r-u_l)*(pl[RHO]+pr[RHO])*(cs_l+cs_r);
  p_s = MAX(p_s, 0.0);

  // if (pratio > 2.0) {
  //   if (p_s < p_min) {
  //     // two rarefaction waves
  //     double Gamma_avg = 0.5*(Gamma_l + Gamma_r);
  //     double plr = pow(p_l/p_r, (Gamma_avg-1.0)/(2.0*Gamma_avg));
  //     ss = (plr*u_l/cs_l + u_r/cs_r + 2.0*(plr - 1.0)/(Gamma_avg - 1.0))/(plr/cs_l + 1.0/cs_r);
  //     p_s = 0.5*(p_l*pow(1.0 + (Gamma_l-1)/(2.0*cs_l)*(u_l - ss),2*Gamma_l/(Gamma_l-1.0)) + p_r*pow(1.0 + (Gamma_r-1.0)/(2.0*cs_r)*(ss - u_r),2.0*Gamma_r/(Gamma_r-1.0)));
  //     p_s = MAX(0.0, p_s);
  //     q_l = (p_s <= p_l) ? 1 : sqrt(1 + (Gamma_l+1)/(2*Gamma_l)*(p_s/p_l - 1.0));
  //     q_r = (p_s <= p_r) ? 1 : sqrt(1 + (Gamma_r+1)/(2*Gamma_r)*(p_s/p_r - 1.0));
  //     sl = u_l - cs_l*q_l;
  //     sr = u_r + cs_r*q_r;
  //   } else {
  //     // two shock waves
  //     double a_l = 2.0/((Gamma_l+1.0)*pl[RHO]);
  //     double a_r = 2.0/((Gamma_r+1.0)*pr[RHO]);
  //     double b_l = (Gamma_l-1.0)/(Gamma_l+1.0)*p_l;
  //     double b_r = (Gamma_r-1.0)/(Gamma_r+1.0)*p_r;
  //     double g_l = sqrt(a_l/(p_s + b_l));
  //     double g_r = sqrt(a_r/(p_s + b_r));
  //     p_s = (g_l*p_l + g_r*p_r - (u_r-u_l))/(g_l + g_r);
  //     p_s = MAX(0.0, p_s);
  //     q_l = (p_s <= p_l) ? 1 : sqrt(1 + (Gamma_l+1)/(2*Gamma_l)*(p_s/p_l - 1.0));
  //     q_r = (p_s <= p_r) ? 1 : sqrt(1 + (Gamma_r+1)/(2*Gamma_r)*(p_s/p_r - 1.0));
  //     sl = u_l - cs_l*q_l;
  //     sr = u_r + cs_r*q_r;
  //   }
  // } else {
  //   q_l = (p_s <= p_l) ? 1 : sqrt(1 + (Gamma_l+1)/(2*Gamma_l)*(p_s/p_l - 1.0));
  //   q_r = (p_s <= p_r) ? 1 : sqrt(1 + (Gamma_r+1)/(2*Gamma_r)*(p_s/p_r - 1.0));
  //   sl = u_l - cs_l*q_l;
  //   sr = u_r + cs_r*q_r;
  // }

  if (pratio < 2.0 && p_s >= p_min && p_s <= p_max) {
    q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
    q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
    // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
    // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
    sl = u_l - cs_l*q_l;
    sr = u_r + cs_r*q_r;
    ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
  } else {
    if (p_s < p_min) {
      // Two rarefaction waves (see Toro, eqs. 9.35-9.36)
      Gamma_avg = 0.5*(Gamma_l + Gamma_r);
      plr = pow(p_l/p_r, (Gamma_avg-1.0)/(2.0*Gamma_avg));
      ss = (plr*u_l/cs_l + u_r/cs_r + 2.0*(plr - 1.0)/(Gamma_avg - 1.0))/(plr/cs_l + 1.0/cs_r);
      p_s = 0.5*(p_l*pow(1.0 + (Gamma_l-1.0)/(2.0*cs_l)*(u_l - ss),2.0*Gamma_l/(Gamma_l-1.0)) + p_r*pow(1.0 + (Gamma_r-1.0)/(2.0*cs_r)*(ss - u_r),2.0*Gamma_r/(Gamma_r-1.0)));
      p_s = MAX(p_s, 0.0);
      q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
      q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
      // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
      // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
      sl = u_l - cs_l*q_l;
      sr = u_r + cs_r*q_r;
      ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
    } else {
      // Two shock waves (see Toro, eqs. 9.42-9.43)
      a_l = 2.0/((Gamma_l+1.0)*pl[RHO]);
      a_r = 2.0/((Gamma_r+1.0)*pr[RHO]);
      b_l = (Gamma_l-1.0)/(Gamma_l+1.0)*p_l;
      b_r = (Gamma_r-1.0)/(Gamma_r+1.0)*p_r;
      g_l = sqrt(a_l/(p_s + b_l));
      g_r = sqrt(a_r/(p_s + b_r));
      p_s = (g_l*p_l + g_r*p_r - (u_r-u_l))/(g_l + g_r);
      p_s = MAX(p_s, 0.0);
      q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
      q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
      // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
      // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
      sl = u_l - cs_l*q_l;
      sr = u_r + cs_r*q_r;
      ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
    }
  }
  
  rhowl = sqrt(pl[RHO]);
  rhowr = sqrt(pr[RHO]);
  SLOOP {
    vriemann[dd] = (rhowl*pl[U1+dd] + rhowr*pr[U1+dd])/(rhowl+rhowr);
  }

  p_uflux(nhydro,pl,g,p_l,dir,ul,fl);
  p_uflux(nhydro,pr,g,p_r,dir,ur,fr);

  cmax = MAX(sr/g->scale[1+dir][dir], 0.0);
  cmin = fabs(MIN(sl/g->scale[1+dir][dir], 0.0));
  *sig_speed = MAX(cmax,cmin);

  HLOOP {
    flux[vv] = area*(cmin*fr[vv] + cmax*fl[vv] - cmax*cmin*(ur[vv]-ul[vv]))/(cmax+cmin);
  }

  TIMER_STOP;

  return 0;
}

int riemann_flux_HLLC(const int nhydro, const int nvars,const double* restrict pl, const double* restrict pr, const zone_geom* restrict g, const int dir, 
                      double* restrict flux, double* restrict sig_speed, double* restrict vriemann)
{
  int vv,dd,ssflag;
  double u_l,u_r;
  double p_l,p_r,cs_l,cs_r,cmax,plr;
  double Gamma_l,Gamma_r;
  double temp_l,temp_r;
  double ent_l,ent_r;
  double sl,sr,ss,p_s,slr;
  double q_l,q_r,pratio,p_max,p_min;
  double a_l,a_r,b_l,b_r,g_l,g_r,z;
  double vsq,Gamma_max,Gamma_min,Gamma_avg;
  double u[nhydro];
  double fx[nvars];
  double ds[nvars];
  double Ul[nvars],Ur[nvars],fl[nvars],fr[nvars];
  double Us[nhydro];
  double vcov[3];
  double area = g->area[dir];
  double eos[NEOS];
  double uint,pnew;
  double ustar;
  double sonicl,sonicr;
  int isshock = 0;

  //TIMER_START("riemann_flux_HLLC");

  vv = nvars;
  p_l = pl[vv];  p_r = pr[vv];  vv++;
  Gamma_l = pl[vv];  Gamma_r = pr[vv];

  if (area < 1.0e-16) {
    *sig_speed = 1.0e-200;
    SLOOP { vriemann[dd] = 0.5*(pl[U1+dd] + pr[U1+dd]); }
    HLOOP { flux[vv] = 0.0; }
    //TIMER_STOP;
    return isshock; 
  }

  u_l = pl[U1+dir]*g->scale[1+dir][dir];
  u_r = pr[U1+dir]*g->scale[1+dir][dir];

  cs_l = sqrt(Gamma_l*p_l/pl[RHO]);
  cs_r = sqrt(Gamma_r*p_r/pr[RHO]);

  if (p_l > p_r) {
    p_max = p_l;
    p_min = p_r;
  } else {
    p_max = p_r;
    p_min = p_l;
  }
  pratio = p_max/p_min;

  isshock = ((p_max-p_min)/p_min > 0.5) ? 1 : 0;

  p_s = 0.5*(p_l + p_r) - 0.125*(u_r-u_l)*(pl[RHO]+pr[RHO])*(cs_l+cs_r);
  p_s = MAX(0,p_s);

  if (pratio < 2.0 && p_s >= p_min && p_s <= p_max) {
    q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
    q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
    // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
    // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
    sl = u_l - cs_l*q_l;
    sr = u_r + cs_r*q_r;
    ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
  } else {
    if (p_s < p_min) {
      // Two rarefaction waves (see Toro, eqs. 9.35-9.36)
      Gamma_avg = 0.5*(Gamma_l + Gamma_r);
      plr = pow(p_l/p_r, (Gamma_avg-1.0)/(2.0*Gamma_avg));
      ss = (plr*u_l/cs_l + u_r/cs_r + 2.0*(plr - 1.0)/(Gamma_avg - 1.0))/(plr/cs_l + 1.0/cs_r);
      p_s = 0.5*(p_l*pow(1.0 + (Gamma_l-1.0)/(2.0*cs_l)*(u_l - ss),2.0*Gamma_l/(Gamma_l-1.0)) + p_r*pow(1.0 + (Gamma_r-1.0)/(2.0*cs_r)*(ss - u_r),2.0*Gamma_r/(Gamma_r-1.0)));
      p_s = MAX(p_s, 0.0);
      q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
      q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
      // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
      // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
      sl = u_l - cs_l*q_l;
      sr = u_r + cs_r*q_r;
      ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
    } else {
      // Two shock waves (see Toro, eqs. 9.42-9.43)
      a_l = 2.0/((Gamma_l+1.0)*pl[RHO]);
      a_r = 2.0/((Gamma_r+1.0)*pr[RHO]);
      b_l = (Gamma_l-1.0)/(Gamma_l+1.0)*p_l;
      b_r = (Gamma_r-1.0)/(Gamma_r+1.0)*p_r;
      g_l = sqrt(a_l/(p_s + b_l));
      g_r = sqrt(a_r/(p_s + b_r));
      p_s = (g_l*p_l + g_r*p_r - (u_r-u_l))/(g_l + g_r);
      p_s = MAX(p_s, 0.0);
      q_l = (p_s <= p_l) ? 1.0 : sqrt(1.0 + (Gamma_l+1.0)/(2.0*Gamma_l)*(p_s/p_l - 1.0));
      q_r = (p_s <= p_r) ? 1.0 : sqrt(1.0 + (Gamma_r+1.0)/(2.0*Gamma_r)*(p_s/p_r - 1.0));
      // sl = MIN(u_l - cs_l*q_l, u_r - cs_r*q_r);
      // sr = MAX(u_l + cs_l*q_l, u_r + cs_r*q_r);
      sl = u_l - cs_l*q_l;
      sr = u_r + cs_r*q_r;
      ss = (p_r - p_l + pl[RHO]*u_l*(sl - u_l) - pr[RHO]*u_r*(sr - u_r))/(pl[RHO]*(sl - u_l) - pr[RHO]*(sr - u_r));
    }
  }

  if (ss > sr) sr = 1.001*ss;
  if (ss < sl) sl = 1.001*ss;

  *sig_speed = MAX(fabs(sl),fabs(sr))/g->scale[1+dir][dir];

  // supersonic cases
  if (sl > 0.0) {
    pflux(nhydro,pl, p_l, g, dir, flux);
    SLOOP {
      vriemann[dd] = pl[U1+dd];
    }
    //TIMER_STOP;
    return isshock;
  } else if (sr < 0.0) {
    pflux(nhydro,pr, p_r, g, dir, flux);
    SLOOP {
      vriemann[dd] = pr[U1+dd];
    }
    //TIMER_STOP;
    return isshock;
  }

  if (ss > 0.0) {
    p_uflux(nhydro,pl, g, p_l, dir, u, fx);
    Us[0] = pl[RHO]*(sl - u_l)/(sl - ss);
    Us[ETOT] = (u[ETOT]*(sl - u_l) - p_l*u_l + p_s*ss)/(sl - ss);
    SLOOP {
      Us[U1+dd] = pl[U1+dd];
    }
    Us[U1+dir] = ss/g->scale[1+dir][dir];
    geom_lower(&Us[U1], vcov, g->gcov[1+dir]);
    vsq = 0.0;
    SLOOP {
      vsq += Us[U1+dd]*vcov[dd];
      vriemann[dd] = Us[U1+dd];
      Us[U1+dd] = vcov[dd]*Us[0];
    }
    for (vv=U1+SPACEDIM;vv<nhydro;vv++) {
      Us[vv] = Us[0]*pl[vv];
    }
    slr = sl/g->scale[1+dir][dir];
  } else {
    p_uflux(nhydro,pr, g, p_r, dir, u, fx);
    Us[0] = pr[RHO]*(sr - u_r)/(sr - ss);
    Us[ETOT] = (u[ETOT]*(u_r - sr) + p_r*u_r - p_s*ss)/(ss - sr);
    SLOOP {
      Us[U1+dd] = pr[U1+dd];
    }
    Us[U1+dir] = ss/g->scale[1+dir][dir];
    geom_lower(&Us[U1], vcov, g->gcov[1+dir]);
    vsq = 0.0;
    SLOOP {
      vsq += Us[U1+dd]*vcov[dd];
      vriemann[dd] = Us[U1+dd];
      Us[U1+dd] = vcov[dd]*Us[0];
    }
    for (vv=U1+SPACEDIM;vv<nhydro;vv++) {
      Us[vv] = Us[0]*pr[vv];
    }
    slr = sr/g->scale[1+dir][dir];
  }

  HLOOP {
    flux[vv] = area*(fx[vv] + slr*(Us[vv] - u[vv]));
  }
  
  //TIMER_STOP;
  return isshock;
}


*/
