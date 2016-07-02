#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F90 = ifort
#F90 = gfortran
#compiling object file flags
OPT1 = -O3 -parallel -openmp -mkl
OPT2 = 
#OPT1 = -O0 -g -CB -warn
#OPT2 = 
FFLAGS = -c $(OPT1) $(OPT2)
#FFLAGS = -c -O0 -g -CB -warn $(OPT2)
#linking flags
LFLAGS = $(OPT1) $(OPT2)
#LFLAGS = -O0 -g -CB -warn $(OPT2)
#testing flags
#LFLAGS = -O0 -g -CB 
#fitsio libraries
FITSIODIR = /usr/local/lib
#Pgplot plot libraries
PGPLOTDIR = /usr/local/lib
#X11 libraries
X11DIR = /usr/X11R6/lib
# libraries for linking PGPLOT
LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng
#LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng 
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/bin/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: psffit makeapmask

apphotincl = precision.o getfits.o readstars.o fluxcoos.o makeapmap.o photometry.o
apphot: apphot.f90 $(apphotincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ apphot.f90 $(apphotincl) $(LIBS2)

psffitincl = precision.o fittermodv3.o getfits.o readstars.o fitterv3.o minpack.o cstarmap.o findnei.o minpack.o fluxcoos.o exportsol.o
psffit: psffit.f90 $(psffitincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ psffit.f90 $(psffitincl) $(LIBS2)

#building object libraries
precision.o: $(UTILS)precision.f90 
	$(F90) $(FFLAGS) $(UTILS)precision.f90
fittermod.o: $(UTILS)fittermod.f90
	$(F90) $(FFLAGS) $(UTILS)fittermod.f90
fittermodv3.o: $(UTILS)fittermodv3.f90
	$(F90) $(FFLAGS) $(UTILS)fittermodv3.f90
getfits.o: $(UTILS)getfits.f90
	$(F90) $(FFLAGS) $(UTILS)getfits.f90
readstars.o: $(UTILS)readstars.f90
	$(F90) $(FFLAGS) $(UTILS)readstars.f90
fitter.o: $(UTILS)fitter.f90
	$(F90) $(FFLAGS) $(UTILS)fitter.f90
fitterv2.o: $(UTILS)fitterv2.f90
	$(F90) $(FFLAGS) $(UTILS)fitterv2.f90
fitterv3.o: $(UTILS)fitterv3.f90
	$(F90) $(FFLAGS) $(UTILS)fitterv3.f90
minpack.o: $(UTILS)minpack.f
	$(F90) $(FFLAGS) $(UTILS)minpack.f
cstarmap.o: $(UTILS)cstarmap.f90
	$(F90) $(FFLAGS) $(UTILS)cstarmap.f90
findnei.o: $(UTILS)findnei.f90
	$(F90) $(FFLAGS) $(UTILS)findnei.f90
lbfgsb.o: $(UTILS)lbfgsb.f
	$(F90) $(FFLAGS) $(UTILS)lbfgsb.f
timer.o: $(UTILS)timer.f
	$(F90) $(FFLAGS) $(UTILS)timer.f
linpack.o: $(UTILS)linpack.f
	$(F90) $(FFLAGS) $(UTILS)linpack.f
fluxcoos.o: $(UTILS)fluxcoos.f90
	$(F90) $(FFLAGS) $(UTILS)fluxcoos.f90
exportsol.o: $(UTILS)exportsol.f90
	$(F90) $(FFLAGS) $(UTILS)exportsol.f90
makeapmap.o: $(UTILS)makeapmap.f90
	$(F90) $(FFLAGS) $(UTILS)makeapmap.f90
photometry.o: $(UTILS)photometry.f90
	$(F90) $(FFLAGS) $(UTILS)photometry.f90

# Removing object files
.PHONY : clean
clean :
	rm *.o
