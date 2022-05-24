# MAKEFILE OPTIONS FOR:  tiger, Intel compilers

# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = h5pcc
F90     = h5pfc
CFLAGS  = -std=c99 -O3 -diag-disable 161,3180,6843
FFLAGS  = -cpp -g -O3 -assume byterecl -diag-disable 3180,6843
LDFLAGS = -lgsl -mkl -lm
