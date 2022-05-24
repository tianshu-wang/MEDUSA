# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC      = h5pcc
F90     = h5pfc
CFLAGS  = -std=c99 -O2 -I/usr/local/include 
#CFLAGS  = -std=c99 -g -Og -Wall -I/usr/local/include
#FFLAGS  = -cpp -Og -g
