# -*- mode: makefile -*-
# -----------------------------------------------------------------
# Makefile for Rooster serial examples
#
# This file can be used as a template for
# other user Makefiles.
# -----------------------------------------------------------------

SHELL = sh

# ATTENTION: path to SUNDIALS installation directory
prefix     = /usr/local
includedir = ${prefix}/fortran
libdir     = ${prefix}/lib

# ATTENTION: path to FORTRAN compiler
F90      = /usr/bin/gfortran
F90FLAGS = -Wall -O2 -fopenmp
F90LIBS  = -lm -fopenmp

# ATTENTION: build target
TARGET   = Rooster_MF

# ------------------------------------------------------------------------------

INCLUDES  = -I${includedir}
LIBRARIES = -lsundials_fcvode_mod -lsundials_cvode -lsundials_fcore_mod ${F90LIBS}
LINKFLAGS = -Wl,-rpath,$(libdir)

# ------------------------------------------------------------------------------

FILES    = A0_Rooster_fmod
FILES_DEPENDENCIES = C2_Glodata_mod \
           C1_Interp_mod C1_Getdata_mod C1_Matrix_mod \
           C0_Correlate_mod \
           B3_Fuel_mod B3_Gasgap_mod B3_Clad_mod B3_Thermbc_mod\
           B2_Fuelrod_mod B2_Htstr_mod B2_Pipe_mod B2_Junction_mod \
           B1_Solid_mod B1_Fluid_mod B1_Signal_mod B1_Core_mod\
           B0_Reactor_mod \
           A2_Coupling_mod A2_Control_mod \
           A1_Comprhs_mod  \

OBJECTS =  ${FILES:=.o}
OBJECTS_DEPENDENCIES = ${FILES_DEPENDENCIES:=.o}

# ------------------------------------------------------------------------------

.SUFFIXES : .o .f90
.f90.o :
	${F90} ${F90FLAGS} ${INCLUDES} -c $<

# ------------------------------------------------------------------------------

all: ${OBJECTS_DEPENDENCIES} ${OBJECTS}
	@for i in ${FILES} ; do \
	  echo "${F90} -o ${TARGET} $${i}.o ${OBJECTS_DEPENDENCIES} ${F90FLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${F90} -o ${TARGET} $${i}.o ${OBJECTS_DEPENDENCIES} ${F90FLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

clean:
	rm -f *.o *.mod
	rm -f ${OBJECTS}
	rm -f ${OBJECTS_DEPENDENCIES}
	rm -f ${FILES}
	rm -f ${FILES_DEPENDENCIES}
# ------------------------------------------------------------------------------
