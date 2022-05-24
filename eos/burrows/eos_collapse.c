
#include "../../decs.h"

#if USE_FORTRAN_CODE
void FORT_EOS_CALL(eos_given_rex)(double *Gamma, double *press, double *cs, double *temp, double *ent, double *dpdr, double *dpde, double *rho, double *e, double *ye, int *pt);
void FORT_EOS_CALL(eos_given_rpx)(double *rho, double *p, double *temp, double *e, double *ye, double *Gamma, int *pt);
void FORT_EOS_CALL(eos_given_rtx)(double *e, double *press, double *gam, double *cs, double *rho, double *temp, double *ye, int *pt);
void FORT_EOS_CALL(eos_given_rsx)(double *gam, double *press, double *cs, double *temp, double *e, double *rho, double *ent, double *ye, int *pt);
void FORT_EOS_CALL(min_e_given_rx)(double *rho, double *ye, double *e, int *pt);
void FORT_EOS_CALL(eos_get_etaxnxp)(double *eta, double *xn, double *xp, double *rho, double *T, double *ye, int *pt);
#else
#include "./eos_stuff.h"
#endif

void eos_given_temp(const double *table,double *p, double *eos)
{
  int *pt_index=NULL;
  double eint;
#if USE_FORTRAN_CODE
  FORT_EOS_CALL(eos_given_rtx)(&eint, &eos[PRESS], &eos[GAMMA], &eos[CS], &p[RHO], &eos[TEMP], &p[YE], pt_index);
#else
  eos_given_rtx(table,&eint, &eos[PRESS], &eos[GAMMA], &eos[CS], &p[RHO], &eos[TEMP], &p[YE]);
#endif

  return;
}

#if(USE_LINEAR_ALLOCATION==TRUE)
int eos_fill_tianshu(const double *table,double **p, double **eos)
#else
int eos_fill_tianshu(const double *table,double NDP_PTR p, double NDP_PTR eos)
#endif
{
  int *pt_index=NULL;
  double Gamma, dpdr, dpde, e;
  double press,cs,tp,s;


  eos_given_rex_tianshu(table,p,eos);


  return 1;
}

#if(USE_LINEAR_ALLOCATION==TRUE)
int eos_fill_ghost_tianshu(const double *table,double **p, double **eos)
#else
int eos_fill_ghost_tianshu(const double *table,double NDP_PTR p, double NDP_PTR eos)
#endif
{
  int *pt_index=NULL;
  double Gamma, dpdr, dpde, e;
  double press,cs,tp,s;


  eos_given_rex_ghost_tianshu(table,p,eos);


  return 1;
}

int eos_fill(const double *table,double *p, double *eos)
{
  int *pt_index=NULL;
  double Gamma, dpdr, dpde, e;
  double press,cs,tp,s;

  e = p[UU]/p[RHO];
  if (eos[TEMP] < 1.0e-10) tp = 1.0;
  else tp = eos[TEMP];
#if USE_FORTRAN_CODE
  FORT_EOS_CALL(eos_given_rex)(&Gamma, &press, &cs, &tp, &s, &dpdr, &dpde, &p[RHO], &e, &p[YE], pt_index);
#else
  eos_given_rex(table,&Gamma, &press, &cs, &tp, &s, &dpdr, &dpde, &p[RHO], &e, &p[YE]);
#endif

  eos[PRESS] = press;
  eos[CS] = cs;
  eos[TEMP] = tp;
  eos[ENT] = s;
  eos[GAMMA] = Gamma;

  return 1;
}

double eos_zero_point(const double *table,double rho, double ye)
{
  int *pt_index=NULL;
  double emin;
#if USE_FORTRAN_CODE
  FORT_EOS_CALL(min_e_given_rx)(&rho, &ye, &emin, pt_index);
#else
  min_e_given_rx(table,rho, ye, &emin);
#endif
  return rho*emin;
}

double eos_rpx(const double *table,double rho, double press, double ye, double Tguess, double *Gamma)
{
  int *pt_index=NULL;
  double e;

  //fprintf(stderr,"calling fort_rpx with: %g %g %g %g\n", rho, press, ye, Tguess);
  //fflush(stderr);
#if USE_FORTRAN_CODE
  FORT_EOS_CALL(eos_given_rpx)(&rho, &press, &Tguess, &e, &ye, Gamma, pt_index);
#else
  eos_given_rpx(table,&rho, &press, &Tguess, &e, &ye, Gamma);
#endif
  //fprintf(stderr,"Energy3 is %g\n", rho*e);

  return rho*e;
}

double eos_get_Temp(const double *table,double rho, double e, double ye, double Tguess)
{
  int *pt_index = NULL;
  double Gamma, press, cs, s, dpdr, dpde;
#if USE_FORTRAN_CODE
  FORT_EOS_CALL(eos_given_rex)(&Gamma, &press, &cs, &Tguess, &s, &dpdr, &dpde, &rho, &e, &ye, pt_index);
#else
  eos_given_rex(table,&Gamma, &press, &cs, &Tguess, &s, &dpdr, &dpde, &rho, &e, &ye);
#endif
  return Tguess;
}

void eos_p_cs(const double *table,double *p, double *press, double *cs, double *Gamma)
{
  int ie;
  double eos[NEOS];

  for (ie=0; ie<NEOS; ie++) eos[ie] = 0.0;
  eos[TEMP] = temp_guess;

  //*Gamma =
  eos_fill(table,p, eos);

  *press = eos[PRESS];
  *cs = eos[CS];
  *Gamma = eos[GAMMA];

  return;
}

void eos_p(const double *table,double *p, double *press)
{
  double cs,Gamma;
  eos_p_cs(table,p, press, &cs, &Gamma);
  return;
}

double eos_cs(const double *table,double *p)
{
  double press,cs,Gamma;
  eos_p_cs(table,p, &press, &cs, &Gamma);
  return cs;
}


double u_given_rpx(const double *table,double rho, double press, double ye, double Tguess)
{
  int *pt_index=NULL;
  double e,Gamma;
#if USE_FORTRAN_CODE
  FORT_EOS_CALL(eos_given_rpx)(&rho, &press, &Tguess, &e, &ye, &Gamma, pt_index);
#else
  eos_given_rpx(table,&rho, &press, &Tguess, &e, &ye, &Gamma);
#endif
  return rho*e;
}

double u_given_rtx(const double *table,double rho, double temp, double ye)
{
  int *pt_index=NULL;
  double u,press,Gamma,cs;
#if USE_FORTRAN_CODE
  FORT_EOS_CALL(eos_given_rtx)(&u, &press, &Gamma, &cs, &rho, &temp, &ye, pt_index);
#else
  eos_given_rtx(table,&u, &press, &Gamma, &cs, &rho, &temp, &ye);
#endif
  return rho*u;
}

void get_lr_eos(const double *table,double rho, double temp, double ye, double *u, double *press, double *cs, double *Gamma)
{
  int *pt_index=NULL;
#if USE_FORTRAN_CODE
  FORT_EOS_CALL(eos_given_rtx)(u, press, Gamma, cs, &rho, &temp, &ye, pt_index);
#else
  eos_given_rtx(table,u, press, Gamma, cs, &rho, &temp, &ye);
#endif
  *u *= rho;

  return;
}

void eos_rsx(const double *table,double rho, double ent, double ye, double *u, double *press, double *cs, double *Gamma, double *temp)
{
  int *pt_index=NULL;
#if USE_FORTRAN_CODE
  FORT_EOS_CALL(eos_given_rsx)(Gamma, press, cs, temp, u, &rho, &ent, &ye, pt_index);
#else
  eos_given_rsx(table,Gamma, press, cs, temp, u, &rho, &ent, &ye);
#endif
  *u *= rho;

  return;
}

GPU_PRAGMA(omp declare target)
void eos_get_etaxnxp(const double *table,double rho, double temp, double ye, double *eta, double *xn, double *xp)
{
  int *pt_index=NULL;
#if USE_FORTRAN_CODE
  FORT_EOS_CALL(eos_get_etaxnxp)(eta, xn, xp, &rho, &temp, &ye, pt_index);
#else
  eos_get_etaxnxp_(table,eta, xn, xp, &rho, &temp, &ye);
#endif
  return;
}
GPU_PRAGMA(omp end declare target)
