
#include "../../decs.h"
#include "../../constants.h"

void eos_p_cs(double *p, double *press, double *cs, double *Gamma)
{
	*press = (gam-1.0)*p[UU];
	*cs = sqrt(gam*(*press)/p[RHO]);
	*Gamma = gam;
  
	return;
}

void eos_p(double *p, double *press)
{
	*press = (gam-1.0)*p[UU];

	return;
}

double eos_cs(double *p)
{
	double press,cs,Gamma;

	eos_p_cs(p, &press, &cs, &Gamma);

	return cs;
}

void eos_given_temp(double *p, double *eos)
{
	eos_p_cs(p, &eos[PRESS], &eos[CS], &eos[GAMMA]);
}

double eos_get_Temp(double rho, double e, double cmp, double Tguess)
{
	return e*(MP+ME)*(gam-1.0)/(2.0*KBOLTZ);	// ASSUMES ionized hydrogen
}

double u_given_rtx(double rho, double temp, double x)
{
	return 2.0*KBOLTZ*temp*rho/(MP+ME)/(gam-1.);
}

int eos_fill(double *p, double *eos)
{
	double press,cs,tp,Gamma;

	eos_p_cs(p, &press, &cs, &Gamma);
	tp = press/(2.0*p[RHO]*KBOLTZ)*(MP+ME);

	eos[PRESS] = press;
	eos[CS   ] = cs;
	eos[TEMP ] = tp;
	eos[GAMMA] = Gamma;

  return 1;
}


