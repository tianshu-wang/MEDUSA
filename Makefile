# PLEASE TRY NOT TO ALTER THIS MAKEFILE!
# Instead of replacing machine.mk and setup.mk with specific filenames,
# create a symbolic link using, e.g., "ln -fs machines/hubert.mk machine.mk"
# However, if you do alter the Makefile, please do not commit/push it to the
# repository, since this affects others.  If you do "git add -u", then please
# do "git reset Makefile" afterwards.

# Compilation defaults (Overwrite these in the machine.mk file, not here!)
USE_MPI   = TRUE
USE_OMP   = FALSE
USE_GPU   = FALSE
CC        = gcc
F90       = gfortran
CFLAGS    = -std=c99 -O3 -ffast-math -ftree-vectorize -fpeel-loops
FFLAGS    = -cpp -O3
LDFLAGS   = -lm
OMPFLAGS  = -fopenmp

# Setup defaults (Overwrite these in the setup.mk file, not here!)
USE_FORTRAN_CODE = FALSE
USE_EXP_INTERP = FALSE
USE_LARGER_STEP = FALSE
USE_LINEAR_ALLOCATION = TRUE
INEL_DT_CONTROL = FALSE
RAD_CLOSURE = M1
RCOORD = SINH_MODIFIED
THCOORD = POLYTH

include machine.mk

SRC=main.c array.c mpi.c utils.c step.c geometry.c rad_fluid.c fluid.c riemann.c reconstruct.c bc.c io.c eos.c timer.c quadrupole.c tracer.c
FSRC=
HDR=decs.h defs.h constants.h timer.h enumerators.h build.h
EXE=fornax

include setup.mk
-include testing.mk

######################
#--------------------#
#  DONT TOUCH BELOW  # [unless you know what you're doing :)]
#--------------------#
######################

SRC+= $(PROB)

VPATH+= problems metrics opacities

ifeq ($(USE_MPI), TRUE)
	SRC+= decompose.c
endif

ifeq ($(USE_OMP), TRUE)
	CFLAGS+= $(OMPFLAGS)
	FFLAGS+= $(OMPFLAGS)
	LDFLAGS+= $(OMPFLAGS)
endif

ifeq ($(GEOM), SPHERICAL)
	SRC+= metric_spherical_general.c
endif

ifeq ($(GEOM), CYLINDRICAL)
	SRC+= metric_cyl_general.c
endif

ifeq ($(GEOM), CARTESIAN)
	SRC+= metric_cartesian.c
endif

ifneq ($(GRAV), NO_GRAV)
	SRC+= gravity.c
endif

ifeq ($(DO_RADIATION), TRUE)
	SRC+= radiation.c
#    LDFLAGS+= -llapack -lblas
	ifeq ($(OPAC), BREMSSTRAHLUNG)
		SRC+= opac_emis_brem.c
	endif
	ifeq ($(OPAC), PURE_SCATT)
		SRC+= opac_pure_scatt.c
	endif
	ifeq ($(OPAC), NULL)
		SRC+= opac_null.c
	endif
	ifeq ($(OPAC), COLLAPSE)
		SRC+= opac_emis_neutrino.c
		ifeq ($(USE_FORTRAN_CODE), TRUE)
			FSRC+= opacity_table_module_F90.f90
			LDFLAGS+= -lhdf5_fortran
			ifneq ($(EOS), COLLAPSE)
				LDFLAGS+= -lgfortran
				#LDFLAGS+= -L/Applications/mesasdk/lib -lgfortran -lquadmath
				#LDFLAGS+= -L/usr/local/Cellar/gcc48/4.8.0/gcc/lib -lgfortran -lquadmath
			endif
        ifeq ($(USE_MPI), TRUE)
#            LDFLAGS+= -lmpichf90
        endif
		else
			SRC+= opacity_table_module.c
		endif
	endif
	#LDFLAGS+= -lclapack
	#LDFLAGS+= -lgfortran -lquadmath -L/Users/jdolence/codes/lapack-3.5.0 -llapack -lblas -L/Users/jdolence/codes/xblas-1.0.248 -lxblas
#	LDFLAGS+= -llapack -lblas
endif

ifeq ($(DO_RADIATION), FALSE)
	NEUTRINO=FALSE
	PHOTON=FALSE
endif

ifeq ($(EOS), COLLAPSE)
	#SRC+= eos_stellar_collapse.c
	#VPATH+= eos/stellarcollapse
	VPATH+= eos/burrows
	SRC+= eos_collapse.c
	ifeq ($(USE_FORTRAN_CODE), TRUE)
		MODS+= table_module.mod
		FSRC+= collapse_F90.f90 eos_stuff_F90.f90
		#LDFLAGS+= -L/usr/local/Cellar/gcc48/4.8.0/gcc/lib -lgfortran -lquadmath
		#LDFLAGS+= -L/Applications/mesasdk/lib -lgfortran -lquadmath
	#LDFLAGS+= -lgfortran
    ifeq ($(DO_RADIATION), FALSE)
        ifeq ($(USE_MPI), TRUE)
#            LDFLAGS+= -lmpichf90
        endif
    endif
	else
		SRC+= eos_stuff.c
	endif
endif
ifeq ($(EOS), GAMMA_LAW)
	SRC+= eos_gamma.c
	VPATH+= eos/gamma
endif
ifeq ($(EOS), POLYTROPIC)
	SRC+= eos_poly.c
	VPATH+= eos/polytropic
endif

ifeq ($(PERTURB), DENSITY_RANDOM)
	USE_AUX_SRC=TRUE
endif
ifeq ($(PERTURB), VELOCITY_SPHERICAL_HARMONIC)
  ifneq ($(strip $(GSL_DIR)),)
        CFLAGS  += -I${GSL_DIR}/include
	LDFLAGS += -L${GSL_DIR}/lib
  endif
	LDFLAGS += -lgsl -lgslcblas
endif
ifeq ($(DO_ROTATION), TRUE)
	GRAV=SPHERICAL_MULTIPOLE_GRAV
endif

DEFS=NDIM USE_MPI USE_OMP USE_GPU GEOM EOS GRAV GR_MONOPOLE PN_POTENTIAL USE_EXT_SRC\
     USE_AUX_SRC DO_HYDRO PRINT_OPAC_DATA DO_RADIATION PHOTON NEUTRINO PERTURB\
     DO_ROTATION ENFORCE_FLOORS DENDRITIC_GRID OUTPUT_PRECISION DO_TIMING DEBUG\
     POLAR_AVG USE_FORTRAN_CODE USE_EXP_INTERP RAD_CLOSURE USE_LARGER_STEP\
     USE_LINEAR_ALLOCATION INEL_DT_CONTROL RCOORD THCOORD

define echo_tool
echo $(1) >> build.h;
endef
define echo_tool_new
echo -n $(1) >> build.h;
endef

OBJ_DIR = obj/
FOBJ = $(addprefix $(OBJ_DIR),$(notdir $(FSRC:.f90=.fo)))
COBJ = $(addprefix $(OBJ_DIR),$(notdir $(SRC:.c=.o)))
OBJ = $(FOBJ) $(COBJ)

.PHONY: all clean whipe

# ensure that build.h is made first if we use -j flag
all: $(OBJ_DIR)
	@$(MAKE) --no-print-directory build.h
	@$(MAKE) --no-print-directory $(EXE)
#	codesign -f -s hydro_cert $(EXE)

# For debugging variables in Makefile, e.g. by "make print-OBJ"
print-%  : ; @echo $* = $($*)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(EXE): $(HDR) $(OBJ) $(SRC) $(FSRC)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

table_module.mod: $(OBJ_DIR)collapse_F90.fo

$(OBJ_DIR)%.fo : %.f90
	$(F90) -c $(FFLAGS) -DNDIM=${NDIM} -DUSE_MPI=${USE_MPI} $< -o $@

$(OBJ_DIR)%.o : %.c $(HDR) | build.h
	$(CC) -c $(CFLAGS) $< -o $@

build.h: $(MODS)
	@echo "Creating build.h with compile time definitions"
	@echo "#ifndef BUILD_H_" >> build.h
	@echo "#define BUILD_H_" >> build.h
	@echo "" >> build.h
	@echo "#include \"enumerators.h\"" >> build.h
	@echo "" >> build.h
	@$(foreach def,$(DEFS),$(call echo_tool,"#define $(def) $($(def))"))
	@echo "" >> build.h
	@echo -n "#define BUILD_STRING " >> build.h
	@echo -n "\"This Fornax executable was built with the below options:\\\n" >> build.h
	@$(foreach def,$(DEFS),$(call echo_tool_new,"$(def)=$($(def))\\\n"))
	@echo -n "\"" >> build.h
	@echo "" >> build.h
	@echo "#endif // BUILD_H_" >> build.h
	@echo "" >> build.h

clean:
	rm -rf $(EXE) *.o *.i *.mod *~ build.h $(OBJ_DIR)*

wipe:
	rm -rf $(EXE) *.o *.i *.mod *~ dumps/* log build.h

new:
	$(MAKE) wipe; $(MAKE)
