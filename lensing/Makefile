#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
FPATH = /opt/intel/compilers_and_libraries_2017.2.163/mac/bin/intel64/
F77 = $(FPATH)ifort
F90 = $(FPATH)ifort
#F77 = gfortran
#compiling object file flags
OPT1 = -O3 -parallel -mkl -ipo
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
X11DIR = /usr/X11/lib
# libraries for linking PGPLOT
LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng
#LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/Kepler/lensing/bin/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: pointsource micro micro1

periodmatch: periodmatch.f90
	$(F77) $(LFLAGS) -o $(BIN)$@ periodmatch.f90

fitincl = dataio.o getfitpars.o ttcor.o transitmodel.o keplerian.o mandelagol.o occultquad.o fittransitmodel2.o minpack.o micromodel.o
microfit: microfit.f $(fitincl)
	$(F77) $(LFLAGS) -o $(BIN)$@ microfit.f $(LIBS) $(fitincl)
pointsource: pointsource.f 
	$(F77) $(LFLAGS) -o $(BIN)$@ pointsource.f $(LIBS)
micro: micro.f
	$(F77) $(LFLAGS) -o $(BIN)$@ micro.f $(LIBS)
micro1: micro1.f
	$(F77) $(LFLAGS) -o $(BIN)$@ micro1.f $(LIBS)



#building object libraries

dataio.o: $(UTILS)dataio.f
	$(F77) $(FFLAGS) $(UTILS)dataio.f
getfitpars.o: $(UTILS)getfitpars.f
	$(F77) $(FFLAGS) $(UTILS)getfitpars.f
ttcor.o: $(UTILS)ttcor.f
	$(F77) $(FFLAGS) $(UTILS)ttcor.f
transitmodel.o: $(UTILS)transitmodel.f
	$(F77) $(FFLAGS) $(UTILS)transitmodel.f
keplerian.o: $(UTILS)keplerian.f
	$(F77) $(FFLAGS) $(UTILS)keplerian.f
mandelagol.o: $(UTILS)mandelagol.f
	$(F77) $(FFLAGS) $(UTILS)mandelagol.f
occultquad.o: $(UTILS)occultquad.f
	$(F77) $(FFLAGS) $(UTILS)occultquad.f
fittransitmodel2.o: $(UTILS)fittransitmodel2.f
	$(F77) $(FFLAGS) $(UTILS)fittransitmodel2.f
minpack.o: $(UTILS)minpack.f
	$(F77) $(FFLAGS) $(UTILS)minpack.f
micromodel.o: $(UTILS)micromodel.f
	$(F77) $(FFLAGS) $(UTILS)micromodel.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
