#MAKEFILE OPTIONS FOR:  tiger, Intel compilers

# SPECIFY PARALLELISM
USE_MPI=TRUE
#basically untested, issues with hdf?
#USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = h5pcc
F90     = h5pfc
CFLAGS  = -std=c99 -O3 -xCORE-AVX512 -qopt-zmm-usage=high -debug -diag-disable 161,3180,6843
FFLAGS  = -cpp -g -O3 -xCORE-AVX512 -assume byterecl -diag-disable 3180,6843
LDFLAGS = -lgsl -mkl -lm -lifcore
OMPFLAGS = -qopenmp
