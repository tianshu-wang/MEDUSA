#include "decs.h"

#if (USE_LINEAR_ALLOCATION==TRUE)
void update_eos(const double *table,double ** p, double ** eos)
#else
void update_eos(const double *table,double NDP_PTR p, double NDP_PTR eos)
#endif
{
	int ii=0,jj=0,kk=0,success;
  
  TIMER_START("update_eos");


  success = eos_fill_tianshu(table,p,eos);

  TIMER_STOP;

  return;
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void update_eos_ghost(const double *table,double ** p, double ** eos)
#else
void update_eos_ghost(const double *table,double NDP_PTR p, double NDP_PTR eos)
#endif
{

  TIMER_START("update_eos_ghost");

  eos_fill_ghost_tianshu(table,p,eos);

  TIMER_STOP;
  return;
}
