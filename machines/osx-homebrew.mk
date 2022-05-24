# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = mpicc
F90     = mpifort
CFLAGS  = -Wall -std=c99 -O2 -ftree-vectorize
FFLAGS  = -cpp -g -O2 -fallow-argument-mismatch
LDFLAGS = -Wall -lm -lhdf5 -lhdf5_fortran -lgsl -lgslcblas -llapack -lblas -lgfortran -lmpifort -lmpi -lpmpi
