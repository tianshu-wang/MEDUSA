# MAKEFILE OPTIONS FOR:  ThetaGPU, CPU/GPU hybrid nodes, NVIDIA HPC SDK compilers

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
LDFLAGS  = -lm -lnvf -lmpi_mpifh -L/lus/theta-fs0/software/thetagpu/hdf5/1.12.0/lib -lhdf5 -v
#OMPFLAGS = -mp=gpu -gpu=cc80,managed,lineinfo -Minfo=mp,accel
OMPFLAGS = -mp=gpu -gpu=cc80,lineinfo -Minfo=mp,accel 
#OMPFLAGS = -mp=gpu -gpu=cc80 -Minfo=mp,accel
