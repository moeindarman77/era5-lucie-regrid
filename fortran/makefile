#=============================================================================
# Makefile for ERA5→T30 regridder
#
# Prerequisites:
#   1. bspline-fortran must be built first:
#      cd ../bspline-fortran && mkdir -p build && cd build
#      cmake .. -DCMAKE_Fortran_COMPILER=ftn && make
#
#   2. Module environment on Derecho:
#      module load ncarenv intel ncarcompilers cray-mpich netcdf cmake
#
# Usage:
#   make                   # build main_regrid.exe
#   make clean             # remove build artifacts
#   make bspline           # build the bspline-fortran dependency
#=============================================================================

# Compiler — Cray wrapper on Derecho; override with: make COMPILE2=ifort
COMPILE2 ?= ftn

# bspline-fortran paths (built via CMake in submodule)
BSPLINE_BUILD = ../bspline-fortran/build
BSPLINE_INC   = $(BSPLINE_BUILD)/lib
BSPLINE_LIB   = $(BSPLINE_BUILD)/lib/libbspline-fortran.a

COMOTT2 = -O2 -traceback -heap-arrays 10 -I$(BSPLINE_INC)
NETCDF_LIBS = -lnetcdff

FILES = \
       mod_utilities.o \
       stringtype.o \
       mod_io.o \
       mod_interp2d.o \
       regrid_era.o

.SUFFIXES: .f .f90 .o

%.o: %.mod

.f90.o:
	$(COMPILE2) $(COMOTT2) -c $<

main_regrid.exe: $(FILES)
	$(COMPILE2) $(COMOTT2) $(FILES) $(BSPLINE_LIB) -o main_regrid.exe $(NETCDF_LIBS)

.PHONY: clean bspline

clean:
	rm -f *.o *.mod main_regrid.exe

bspline:
	@echo "Building bspline-fortran..."
	mkdir -p $(BSPLINE_BUILD)
	cd $(BSPLINE_BUILD) && cmake .. -DCMAKE_Fortran_COMPILER=$(COMPILE2) && $(MAKE)

mod_utilities.o        : mod_utilities.f90
stringtype.o           : stringtype.f90
mod_io.o               : mod_io.f90 mod_utilities.o stringtype.o
mod_interp2d.o         : mod_interp2d.f90 mod_utilities.o
regrid_era.o           : regrid_era.f90 mod_utilities.o mod_interp2d.o mod_io.o
