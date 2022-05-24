# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = h5pcc
F90     = h5pfc
CFLAGS  = -std=c99 -O2 -ftree-vectorize -fpeel-loops
FFLAGS  = -cpp -g -O2
LDFLAGS =  -lgsl -lgslcblas -lblas -llapack -lhdf5_fortran -lhdf5 -lmpi_usempi -lmpi_mpifh -lmpi -lgfortran -lm
