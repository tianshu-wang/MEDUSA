# MAKEFILE OPTIONS FOR:  GPU on stellar

# SPECIFY PARALLELISM
#USE_MPI=FALSE # basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=TRUE
#USE_OMP=FALSE
USE_GPU = TRUE

CC       = mpicc
F90      = mpif90
CFLAGS   = -O3 -std=c99 --diag_suppress=177,550 --display_error_number
FFLAGS   = -O3 -cpp
LDFLAGS  = -lm -lnvf -lmpi_mpifh -L/usr/local/hdf5/nvhpc-21.5/openmpi-4.1.1/1.10.6/lib64 -lhdf5 -v
OMPFLAGS = -mp=gpu -gpu=cc80,lineinfo -Minfo=mp,accel 
