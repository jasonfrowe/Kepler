#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F77 = ifort
F90 = ifort
#compiling object file flags
OPT1 = -O3 -fast -parallel -ipo -no-prec-div
OPT2 = 
FFLAGS = -c $(OPT1) $(OPT2)
#FFLAGS = -c -O0 -g -CB $(OPT2)
#linking flags
LFLAGS = $(OPT1) $(OPT2)
#LFLAGS = -O0 -g -CB $(OPT2)
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
BIN = /Users/rowe/Documents/bin/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: transitfind tpsclean tpsstats effplot transitfind2

pickincl = precision.o
pickcandidates: pickcandidates.f90 $(pickincl)
	$(F90) $(LFLAGS) -o /Users/rowe/Documents/Kepler/periodsearch/pickcandidates pickcandidates.f90 $(pickincl)

find2incl = precision.o rqsort.o stdev.o graphics.o
transitfind2: transitfind2.f90 $(find2incl)
	$(F90) $(LFLAGS) -o $(BIN)transitfind2 transitfind2.f90 $(find2incl) $(LIBS)

notesmerge: notesmerge.f
	$(F77) $(LFLAGS) -o /home/rowe/Kepler/tpsMQ/q1q10/notes/work/notesmerge notesmerge.f

tpsdawg: tpsdawg.f
	$(F77) $(LFLAGS) -o /home/rowe/Kepler/tpsMQ/tpsdawg tpsdawg.f $(LIBS)

effplot: effplot.f
	$(F77) $(LFLAGS) -o $(BIN)effplot effplot.f $(LIBS)

tpsstats: tpsstats.f
	$(F77) $(LFLAGS) -o $(BIN)tpsstats tpsstats.f $(LIBS)

tpsclean: tpsclean.f 
	$(F77) $(LFLAGS) -o $(BIN)tpsclean tpsclean.f

FINDINCL = dataio.o bls.o nyquest.o rqsort.o
transitfind: transitfind1.f $(FINDINCL)
	$(F77) $(LFLAGS) -o $(BIN)transitfind transitfind1.f $(FINDINCL) $(LIBS)

#building object libraries

graphics.o: $(UTILS)graphics.f
	$(F77) $(FFLAGS) $(UTILS)graphics.f
stdev.o: $(UTILS)stdev.f
	$(F77) $(FFLAGS) $(UTILS)stdev.f
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
dataio.o: $(UTILS)dataio.f
	$(F77) $(FFLAGS) $(UTILS)dataio.f
bls.o: $(UTILS)bls.f
	$(F77) $(FFLAGS) $(UTILS)bls.f
nyquest.o: $(UTILS)nyquest.f
	$(F77) $(FFLAGS) $(UTILS)nyquest.f
rqsort.o: $(UTILS)rqsort.f
	$(F77) $(FFLAGS) $(UTILS)rqsort.f

# Removing object files
.PHONY : clean
clean :
	rm *.o ; rm *.mod
