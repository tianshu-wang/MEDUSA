# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = h5pcc
F90     = h5pfc
CFLAGS  = -g -std=c99 -O3 -ftree-vectorize -ffast-math -fno-math-errno -fno-signed-zeros -fno-trapping-math -Wall
FFLAGS  = -cpp -O3
LDFLAGS = -llapack -lblas -lhdf5 -lgfortran -lmpi_mpifh
