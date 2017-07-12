#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F77 = ifort
F90 = ifort
#F77 = gfortran
#compiling object file flags
OPT1 = -O3 -parallel -ipo -openmp -mkl
OPT2 = -heap-arrays 0
#OPT1 = -O3 
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
BIN = /Users/rowe/Documents/Kepler/timing/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: timingpap timingdist

timingfigureincl = 
timingfigure: timingfigure.f90 $(timingfigureincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ timingfigure.f90 $(timingfigureincl) $(LIBS)

timingdistincl = precision.o
timingdist: timingdist.f90 $(timingdistincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ timingdist.f90 $(timingdistincl) $(LIBS)


timingpapincl = precision.o readttfile.o chisqtest.o rqsort.o
timingpap: timingpap.f90 $(timingpapincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ timingpap.f90 $(timingpapincl)

#building object libraries
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
readttfile.o: $(UTILS)readttfile.f90
	$(F90) $(FFLAGS) $(UTILS)readttfile.f90
chisqtest.o: $(UTILS)chisqtest.f90
	$(F90) $(FFLAGS) $(UTILS)chisqtest.f90
rqsort.o: $(UTILS)rqsort.f
	$(F90) $(FFLAGS) $(UTILS)rqsort.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
