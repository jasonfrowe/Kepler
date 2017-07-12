#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F90 = ifort
#F90 = gfortran
#compiling object file flags
OPT1 = -O3 -fast -parallel -ipo -no-heap-arrays
#OPT1 = -O0 -g -CB
#OPT1 = -O3 
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
PLPLOTDIR = -I/usr/local/Cellar/plplot/5.9.11/lib/fortran/modules/plplot -I/usr/local/Cellar/plplot/5.9.11/include/plplot -L/usr/local/Cellar/plplot/5.9.11/lib 
LIBS = $(PLPLOTDIR) -lplplotf95d -lplplotf95cd
#LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng 
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/Kepler/Neptune/bin/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: kpixread

DFITINCL = precision.o fittingmod.o inverse.o pixelfix.o normaldataset.o
datafit: datafit.f90 $(DFITINCL)
	$(F90) $(LFLAGS) -mkl -o $(BIN)$@ datafit.f90 $(LIBS2) $(DFITINCL)

KPIXINCL = precision.o fitsread.o readheader.o getnsteps.o displayfits.o rqsort.o stdev.o writefits.o deletefile.o
kpixread: kpixread.f90 $(KPIXINCL)
	$(F90) $(LFLAGS) -o $(BIN)$@ kpixread.f90 $(LIBS2) $(KPIXINCL)


#building object libraries
normaldataset.o: $(UTILS)normaldataset.f90
	$(F90) $(FFLAGS) $(UTILS)normaldataset.f90
mrqmin.o: $(UTILS)mrqmin.f
	$(F90) $(FFLAGS) $(UTILS)mrqmin.f
pixelfix.o: $(UTILS)pixelfix.f90
	$(F90) $(FFLAGS) $(UTILS)pixelfix.f90
inverse.o: $(UTILS)inverse.f90
	$(F90) $(FFLAGS) $(UTILS)inverse.f90
fittingmod.o: $(UTILS)fittingmod.f90
	$(F90) $(FFLAGS) $(UTILS)fittingmod.f90
minpack.o: $(UTILS)minpack.f
	$(F90) $(FFLAGS) $(UTILS)minpack.f
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
fitsread.o: $(UTILS)fitsread.f90
	$(F90) $(FFLAGS) $(UTILS)fitsread.f90
readheader.o: $(UTILS)readheader.f
	$(F90) $(FFLAGS) $(UTILS)readheader.f
getnsteps.o: $(UTILS)getnsteps.f90
	$(F90) $(FFLAGS) $(UTILS)getnsteps.f90
displayfits.o: $(UTILS)displayfits.f90
	$(F90) $(FFLAGS) $(UTILS)displayfits.f90
rqsort.o: $(UTILS)rqsort.f
	$(F90) $(FFLAGS) $(UTILS)rqsort.f
stdev.o: $(UTILS)stdev.f
	$(F90) $(FFLAGS) $(UTILS)stdev.f
writefits.o: $(UTILS)writefits.f90
	$(F90) $(FFLAGS) $(UTILS)writefits.f90
deletefile.o: $(UTILS)deletefile.f
	$(F90) $(FFLAGS) $(UTILS)deletefile.f
	
# Removing object files
.PHONY : clean
clean :
	rm *.o
	rm *.mod
