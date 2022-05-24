# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = h5pcc
F90     = h5pfc
CFLAGS  = -std=c99 -g -O2 -xCORE-AVX2
FFLAGS  = -cpp -g -O2 -xCORE-AVX2 -assume byterecl
LDFLAGS = -lgsl -mkl -lhdf5_fortran -lifcore -lmpifort -lmpi -lmpigi -ldl -lrt -lpthread
