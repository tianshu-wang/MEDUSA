# MAKEFILE OPTIONS FOR:  cori, KNL nodes, Intel compilers

# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE  # DO NOT ENABLE THIS UNTIL FURTHER NOTICE

CC       = cc
F90      = ftn
CFLAGS   = -O3 -std=c99 -xMIC-AVX512 -diag-disable 161,3180,6843
FFLAGS   = -O3 -cpp -xMIC-AVX512 -diag-disable 3180,6843 -assume byterecl
LDFLAGS  = -lm -L/usr/common/software/gsl/2.1/intel/lib -lgsl -lgslcblas
OMPFLAGS = -qopenmp

