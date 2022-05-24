# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = mpicc
F90     = mpifort
CFLAGS  = $(LOCAL_INCLUDE) -Wall -std=c99 -O2 -ftree-vectorize -fpeel-loops
FFLAGS  = $(LOCAL_FFLAGS) -cpp -g -O2
LDFLAGS = $(LOCAL_LDFLAGS) -Wall -lm -llapack -lblas -lhdf5 -lhdf5_fortran -lgsl -lgslcblas -lgfortran -lmpi_mpifh -lmpi_usempif08
