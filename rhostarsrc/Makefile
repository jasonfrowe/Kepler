#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F77 = ifort
F90 = ifort
#compiling object file flags
OPT1 = -O3 -fast -parallel -ipo -no-prec-div
OPT2 = 
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
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/smodels/YY/
BINISO = /Users/rowe/Documents/smodels/dartmouth/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: rhostar rhoboot2 rhoboot rhorand rhoboot3

RHOSTARISOINCL = precision.o isochrones.o
rhostariso: rhostariso.f90 $(RHOSTARISOINCL)
	$(F90) $(LFLAGS) -o $(BINISO)$@ rhostariso.f90 $(RHOSTARISOINCL)

rhoboot3: rhoboot3.f
	$(F77) $(LFLAGS) -o $(BIN)$@ rhoboot3.f $(LIBS)

rhorand: rhorand.f
	$(F77) $(LFLAGS) -o $(BIN)$@ rhorand.f

rhorand314: rhorand314.f
	$(F77) $(LFLAGS) -o /Users/rowe/Documents/transitfit/Kepler_n/koi314.n/$@ rhorand314.f

RHOSTARINCL = getinputpars.o gettrack.o yytrack.o getage.o spline.o random.o mcmc.o getrhosig.o lininterp.o mcoresig.o massscan.o
rhostar: rhostar.f $(RHOSTARINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ rhostar.f $(RHOSTARINCL)

rhoboot: rhoboot.f
	$(F77) $(LFLAGS) -o $(BIN)$@ rhoboot.f $(LIBS)

rhoboot2: rhoboot2.f
	$(F77) $(LFLAGS) -o $(BIN)$@ rhoboot2.f $(LIBS)

FITINCL = gettrack.o yytrack.o minpack.o lininterp.o mrqmin.o
rhostarfit: rhostarfit.f $(FITINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ rhostarfit.f $(FITINCL)



#building object libraries

isochrones.o: $(UTILS)isochrones.f90
	$(F90) $(FFLAGS) $(UTILS)isochrones.f90
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
mrqmin.o: $(UTILS)mrqmin.f
	$(F77) $(FFLAGS) $(UTILS)mrqmin.f
minpack.o: $(UTILS)minpack.f
	$(F77) $(FFLAGS) $(UTILS)minpack.f
massscan.o: $(UTILS)massscan.f
	$(F77) $(FFLAGS) $(UTILS)massscan.f
mcoresig.o: $(UTILS)mcoresig.f
	$(F77) $(FFLAGS) $(UTILS)mcoresig.f
lininterp.o: $(UTILS)lininterp.f
	$(F77) $(FFLAGS) $(UTILS)lininterp.f
getrhosig.o: $(UTILS)getrhosig.f
	$(F77) $(FFLAGS) $(UTILS)getrhosig.f
yytrack.o: $(UTILS)yytrack.f
	$(F77) $(FFLAGS) $(UTILS)yytrack.f
getinputpars.o: $(UTILS)getinputpars.f
	$(F77) $(FFLAGS) $(UTILS)getinputpars.f
gettrack.o: $(UTILS)gettrack.f
	$(F77) $(FFLAGS) $(UTILS)gettrack.f
getage.o: $(UTILS)getage.f
	$(F77) $(FFLAGS) $(UTILS)getage.f
spline.o: $(UTILS)spline.f
	$(F77) $(FFLAGS) $(UTILS)spline.f
random.o: $(UTILS)random.f
	$(F77) $(FFLAGS) $(UTILS)random.f
mcmc.o: $(UTILS)mcmc.f
	$(F77) $(FFLAGS) $(UTILS)mcmc.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
