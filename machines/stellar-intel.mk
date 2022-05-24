# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
USE_GPU=FALSE
#USE_OMP=TRUE

CC      = h5pcc
F90     = h5pfc
CFLAGS  = -std=c99 -g -O2 -xCORE-AVX512 -diag-disable=161,3180
FFLAGS  = -cpp -g -O2 -xCORE-AVX512 -assume byterecl -diag-disable=6843
LDFLAGS = -lgsl -mkl -lhdf5_fortran -lifcore -lmpifort -lmpi -ldl -lrt -lpthread
