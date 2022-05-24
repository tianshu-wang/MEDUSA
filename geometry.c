
#include "decs.h"

void init_geometry()
{
  TIMER_START("init_geometry");

  if (mpi_io_proc()) {
    fprintf(stderr,"Entering init_geometry...");
    fflush(stderr);
  }

  if (mpi_io_proc()) {
    fprintf(stderr,"1...");
    fflush(stderr);
  }

  init_interp_vol();
  if (mpi_io_proc()) {
    fprintf(stderr,"2...");
    fflush(stderr);
  }

  init_volume();
  if (mpi_io_proc()) {
    fprintf(stderr,"3...");
    fflush(stderr);
  }

  init_area();
  if (mpi_io_proc()) {
    fprintf(stderr,"4...");
    fflush(stderr);
  }

  init_gcov();
  if (mpi_io_proc()) {
    fprintf(stderr,"5...");
    fflush(stderr);
  }

  init_gcon();
  if (mpi_io_proc()) {
    fprintf(stderr,"6...");
    fflush(stderr);
  }

  init_scale();
  if (mpi_io_proc()) {
    fprintf(stderr,"7...");
    fflush(stderr);
  }

  init_conn();

  if (mpi_io_proc()) {
    fprintf(stderr,"done!\n");
    fflush(stderr);
  }

  TIMER_STOP;
  return;
}

double geom_dot(const double* restrict v1, const double* restrict v2, const double* restrict gcov)
{
  int i;
  double vsq = 0.0;

  for (i=0; i<SPACEDIM; i++) vsq += v1[i]*v2[i]*gcov[i];

  return vsq;
}

void geom_lower(const double* restrict vcon, double* restrict vcov, const double* restrict  gcov)
{
  int i;  

  for (i=0; i<SPACEDIM; i++)  vcov[i] = gcov[i]*vcon[i];

  return;
}

void geom_raise(const double* restrict vcov, double* restrict vcon, const double* restrict gcon)
{
  int i; 
  
  for (i=0; i<SPACEDIM; i++)  vcon[i] = gcon[i]*vcov[i];

  return;
}
