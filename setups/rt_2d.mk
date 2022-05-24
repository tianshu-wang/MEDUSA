# MAKEFILE OPTIONS FOR:  init_rt, 2-d 
# (Note, must currently use '3d' branch!)

# SPECIFY DEBUGGING FLAGS
DEBUG=FALSE
PRINT_OPAC_DATA=FALSE

# SPECIFY THE NUMBER OF DIMENSIONS
NDIM=2

# SPECIFY THE PROBLEM FILE
PROB=init_rt.c

# SPECIFY THE GEOMETRY AND GRID STRUCTURE
GEOM=CARTESIAN
DENDRITIC_GRID=FALSE

# CHECK ENERGY AND DENSITY FLOORS?
ENFORCE_FLOORS=FALSE

# USE POLAR AVERAGING SCHEME
POLAR_AVG=FALSE

#---------------------#
# SPECIFY THE PHYSICS
#---------------------#
# HYDRO?
DO_HYDRO=TRUE

# GRAVITY?
PN_POTENTIAL=FALSE
GRAV=FIXED_GRAV
GR_MONOPOLE=FALSE

# RADIATION?
DO_RADIATION=FALSE
NEUTRINO=FALSE  # EITHER NEUTRINO OR PHOTON, BUT NOT BOTH
PHOTON=FALSE

# OPACITIES/EMISSIVITIES
OPAC=NONE

# PERTURBATIONS? (See Muller & Janka 2015, Summa et al. 2016)
# Make sure libgsl and libgslcblas can be found, or else set GSL_DIR somewhere
PERTURB=NONE

# ROTATING MODEL?
DO_ROTATION=FALSE

# SHOULD WE CALL USER DEFINED SOURCE FUNCTIONS?
USE_EXT_SRC=FALSE
USE_AUX_SRC=FALSE

# EOS
EOS=GAMMA_LAW

# DUMP SINGLE OR DOUBLE PRECISION?
OUTPUT_PRECISION=SINGLE

# SHOULD WE TIME THE RUN?
DO_TIMING=TRUE