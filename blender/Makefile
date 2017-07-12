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
BIN = /home/rowe/p555/smodels/YY/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: handmixer

HANDMIXERINCL = getfitpars.o gettrack.o lininterp.o yytrack.o likelihood.o transitmodel.o keplerian.o eclmod.o mandelagol.o dataio.o cfunk.o fittransitmodel2.o mrqmin.o limbdarken.o
handmixer: handmixer.f $(HANDMIXERINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ handmixer.f $(HANDMIXERINCL)

#building object libraries

limbdarken.o: $(UTILS)limbdarken.f
	$(F77) $(FFLAGS) $(UTILS)limbdarken.f
mrqmin.o: $(UTILS)mrqmin.f
	$(F77) $(FFLAGS) $(UTILS)mrqmin.f
cfunk.o: $(UTILS)cfunk.f
	$(F77) $(FFLAGS) $(UTILS)cfunk.f
fittransitmodel2.o: $(UTILS)fittransitmodel2.f
	$(F77) $(FFLAGS) $(UTILS)fittransitmodel2.f
dataio.o: $(UTILS)dataio.f
	$(F77) $(FFLAGS) $(UTILS)dataio.f
mandelagol.o: $(UTILS)mandelagol.f
	$(F77) $(FFLAGS) $(UTILS)mandelagol.f
eclmod.o: $(UTILS)eclmod.f
	$(F77) $(FFLAGS) $(UTILS)eclmod.f
keplerian.o: $(UTILS)keplerian.f
	$(F77) $(FFLAGS) $(UTILS)keplerian.f
transitmodel.o: $(UTILS)transitmodel.f
	$(F77) $(FFLAGS) $(UTILS)transitmodel.f
likelihood.o: $(UTILS)likelihood.f
	$(F77) $(FFLAGS) $(UTILS)likelihood.f
gettrack.o: $(UTILS)gettrack.f	
	$(F77) $(FFLAGS) $(UTILS)gettrack.f
getfitpars.o: $(UTILS)getfitpars.f
	$(F77) $(FFLAGS) $(UTILS)getfitpars.f
lininterp.o: $(UTILS)lininterp.f
	$(F77) $(FFLAGS) $(UTILS)lininterp.f
yytrack.o: $(UTILS)yytrack.f
	$(F77) $(FFLAGS) $(UTILS)yytrack.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
