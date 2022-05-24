# MAKEFILE OPTIONS FOR:  cori, KNL nodes, GNU compilers

# SPECIFY PARALLELISM
#USE_MPI=FALSE basically untested, issues with hdf?
USE_MPI=TRUE
USE_OMP=FALSE
#USE_OMP=TRUE

CC       = cc
F90      = ftn
CFLAGS   = -g -std=c99 -O3 -ffast-math -ftree-vectorize -funroll-loops -ftree-vectorizer-verbose=2 -march=knl
FFLAGS   = -g -cpp -O3 -march=knl
LDFLAGS  = -lm
OMPFLAGS = -fopenmp

