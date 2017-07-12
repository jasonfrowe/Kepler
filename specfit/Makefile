#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F90 = ifort
F77 = ifort
#compiling object file flags
OPT1 = -O3 -fast -parallel -ipo -no-prec-div
OPT2 = 
#OPT1 = -O0 -g -CB -warm
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

all: specfit

specfitincl = precision.o fittingmod.o readspec.o plotspec.o readmodel.o plotmodel.o binmodel.o rqsort.o fitspectrum.o minpack.o
specfit: specfit.f90 $(specfitincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ specfit.f90 $(LIBS) $(specfitincl)


#building object libraries
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
readspec.o: $(UTILS)readspec.f90
	$(F90) $(FFLAGS) $(UTILS)readspec.f90
plotspec.o: $(UTILS)plotspec.f90
	$(F90) $(FFLAGS) $(UTILS)plotspec.f90
readmodel.o: $(UTILS)readmodel.f90
	$(F90) $(FFLAGS) $(UTILS)readmodel.f90
plotmodel.o: $(UTILS)plotmodel.f90
	$(F90) $(FFLAGS) $(UTILS)plotmodel.f90
binmodel.o: $(UTILS)binmodel.f90
	$(F90) $(FFLAGS) $(UTILS)binmodel.f90
rqsort.o: $(UTILS)rqsort.f
	$(F77) $(FFLAGS) $(UTILS)rqsort.f
opkdmain.o: $(UTILS)opkdmain.f
	$(F77) $(FFLAGS) $(UTILS)opkdmain.f
opkda1.o: $(UTILS)opkda1.f
	$(F77) $(FFLAGS) $(UTILS)opkda1.f
opkda2.o: $(UTILS)opkda2.f
	$(F77) $(FFLAGS) $(UTILS)opkda2.f
fitspectrum.o: $(UTILS)fitspectrum.f90
	$(F90) $(FFLAGS) $(UTILS)fitspectrum.f90
minpack.o: $(UTILS)minpack.f
	$(F77) $(FFLAGS) $(UTILS)minpack.f
fittingmod.o: $(UTILS)fittingmod.f90
	$(F90) $(FFLAGS) $(UTILS)fittingmod.f90

# Removing object files
.PHONY : clean
clean :
	rm *.o
	rm *.mod
