# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = h5pcc
F90     = h5pfc
CFLAGS  = $(LOCAL_INCLUDE) -O0 -g -Wall -Werror -std=c99
FFLAGS  = $(LOCAL_FFLAGS)  -O0 -g -cpp
LDFLAGS = $(LOCAL_LDFLAGS) -Wall -lm -llapack -lblas -lhdf5 -lhdf5_fortran -lgsl -lgslcblas -lgfortran -lmpi_mpifh
