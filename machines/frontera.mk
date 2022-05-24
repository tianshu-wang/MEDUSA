# MAKEFILE OPTIONS FOR:  stampede, KNL nodes, Intel compilers

# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC       = mpicc
F90      = mpif90
CFLAGS   = -O3 -std=c99 -xCORE-AVX512 -diag-disable 161,3180,6843 -I${TACC_HDF5_INC} -I${TACC_GSL_INC}
FFLAGS   = -O3 -fpp -xCORE-AVX512 -diag-disable 3180,6843 -assume byterecl -I${TACC_HDF5_INC}
LDFLAGS  = -lm -L${TACC_HDF5_LIB} -lhdf5 -L${TACC_GSL_LIB} -lgsl -lifcore -mkl
OMPFLAGS = -qopenmp
