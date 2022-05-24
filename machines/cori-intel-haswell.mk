# MAKEFILE OPTIONS FOR:  cori, Haswell nodes, Intel compilers

# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC       = cc
F90      = ftn
CFLAGS   = -g -std=c99 -diag-disable 161,3180,6843
FFLAGS   = -g -cpp -diag-disable 3180,6843 -assume byterecl
LDFLAGS  = -lm 
OMPFLAGS = -qopenmp
