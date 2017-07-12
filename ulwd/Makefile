#Name of C compiler
CCMPL = icc
#Name of Fortran compiler
F77 = ifort
#compiling object file flags
OPT1 = -O3 -ftz -ip -ipo -march=core2 -mtune=core2 -assume buffered_io -warn
OPT2 = -warn alignments -align all -align rec8byte -warn
FFLAGS = -c $(OPT1) $(OPT2)
#FFLAGS = -c -O0 -g -CB -warn $(OPT2)
#linking flags
LFLAGS = $(OPT1) $(OPT2)
#LFLAGS = -O0 -g -CB -warn $(OPT2)
#testing flags
#LFLAGS = -O0 -g -CB 
#fitsio libraries
FITSIODIR = /usr/local/fitsio/lib
#Pgplot plot libraries
PGPLOTDIR = /usr/local/pgplot
#X11 libraries
X11DIR = /usr/X11R6/lib
# libraries for linking PGPLOT
LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng -lsvml
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng -lsvml
#Directory where executable are placed
BIN = /home/rowe/p555/whitedwarf/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: massrad wdchandra

wdchandra: wdchandra.f rkqs.o
	$(F77) $(LFLAGS) -o $(BIN)$@ wdchandra.f rkqs.o

massrad: massrad.f
	$(F77) $(LFLAGS) -o $(BIN)$@ massrad.f $(LIBS)

#building object libraries

rkqs.o: $(UTILS)rkqs.f
	$(F77) $(FFLAGS) $(UTILS)rkqs.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
