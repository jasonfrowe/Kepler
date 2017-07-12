#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F77 = ifort
F90 = ifort
#F77 = gfortran
#compiling object file flags
OPT1 = -O3 -fast -parallel -ipo -no-prec-div
OPT2 = 
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
PLPLOTDIR = -I/usr/local/Cellar/plplot/5.9.11/lib/fortran/modules/plplot -I/usr/local/Cellar/plplot/5.9.11/include/plplot -L/usr/local/Cellar/plplot/5.9.11/lib 
#LIBS = $(PLPLOTDIR) -lplplotf95d -lplplotf95cd
LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng 
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/Kepler/Q1Q12/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: catplot

catplotincl = 
catplot: catplot.f90 $(catplotincl)
	$(F90) $(LFLAGS) -o $(BIN)bin/$@ catplot.f90 $(LIBS) $(catplotincl)

q1q12nlistincl = 
q1q12nlist: $(q1q12nlistincl)
	$(F77) $(LFLAGS) -o $(BIN)models/$@ q1q12nlist.f $(q1q12nlistincl)

#building object libraries
histogram2.o: $(UTILS)histogram2.f
	$(F77) $(FFLAGS) $(UTILS)histogram2.f
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90

# Removing object files
.PHONY : clean
clean :
	rm *.o
