# MAKEFILE OPTIONS FOR:  bethe, Intel compilers
# Note that the current modules use a very old version of GCC (4.4.7 from 2012)
# which does not contain necessary language features from OpenMP v3.1 on, so
# do not attempt to compile with USE_OMP=TRUE

# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = h5pcc
F90     = h5pfc
CFLAGS  = -std=c99 -O3 -ftree-vectorize -fpeel-loops
FFLAGS  = -cpp -g -O3
LDFLAGS = -lm -llapack -lblas -lhdf5 -lhdf5_fortran -lgsl -lgslcblas
