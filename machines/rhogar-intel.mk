# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = /opt/intel/compilers_and_libraries_2020.1.217/linux/bin/intel64/icc
F90     = h5pfc
CFLAGS  = $(LOCAL_INCLUDE) -std=c99 -g -O2 -xCORE-AVX512 -diag-disable=161,3180
FFLAGS  = $(LOCAL_FFLAGS) -cpp -g -O2 -xCORE-AVX512 -assume byterecl -diag-disable=6843
LDFLAGS = $(LOCAL_LDFLAGS) -lgsl -mkl -lifcore -lmpi -ldl -lrt -lpthread -lhdf5
