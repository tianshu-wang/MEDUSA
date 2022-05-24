# MAKEFILE OPTIONS FOR:  tiger, GNU compilers

# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = h5pcc
F90     = h5pfc
CFLAGS  = -std=c99 -O3 -ftree-vectorize -fpeel-loops
FFLAGS  = -cpp -g -O3
LDFLAGS = -L/home/mskinner/TOOLS/gsl-2.1 -lgsl -lgslcblas -lblas -llapack -lhdf5_fortran -lhdf5
