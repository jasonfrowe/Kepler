#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F77 = ifort
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
FITSIODIR = /usr/local/fitsio/lib
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

all: transitfit4 transitmcmc mcmchist stellardensity p2g transitchisq4

timingplots: timingplots.f
	$(F77) $(LFLAGS) -o $(BIN)$@ timingplots.f $(LIBS)
chisqincl = dataio.o getfitpars.o transitmodel.o mandelagol.o keplerian.o occultquad.o
transitchisq4: transitchisq4.f $(chisqincl)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitchisq4.f $(chisqincl)

stellardensity: stellardensity.f getfitpars.o
	$(F77) $(LFLAGS) -o $(BIN)$@ stellardensity.f getfitpars.o 

P2GINCL = dataio.o
p2g: p2g.f $(P2GINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ p2g.f $(P2GINCL) 

mcmchist: mcmchist.f 
	$(F77) $(LFLAGS) -o $(BIN)$@ mcmchist.f $(LIBS)

MCMCINCL = getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o mcmc.o ran2.o occultquad.o
transitmcmc: transitmcmc.f $(MCMCINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitmcmc.f $(MCMCINCL)

TRANSINCL3 = getfitpars.o transitmodel.o mandelagol.o dataio.o fittransitmodel.o mrqmin.o keplerian.o lininterp.o readttfile.o occultquad.o
transitfit4: transitfit4.f $(TRANSINCL3)
	$(F77) $(LFLAGS) -o $(BIN)$@ \
	transitfit4.f $(TRANSINCL3) $(LIBS)
transittiming: transittiming.f $(TRANSINCL3)
	$(F77) $(LFLAGS) -o $(BIN)$@ transittiming.f $(TRANSINCL3) $(LIBS)


#building object libraries

occultquad.o: $(UTILS)occultquad.f
	$(F77) $(FFLAGS) $(UTILS)occultquad.f
readttfile.o: $(UTILS)readttfile.f
	$(F77) $(FFLAGS) $(UTILS)readttfile.f
lininterp.o: $(UTILS)lininterp.f
	$(F77) $(FFLAGS) $(UTILS)lininterp.f
getfitpars.o: $(UTILS)getfitpars.f
	$(F77) $(FFLAGS) $(UTILS)getfitpars.f
transitmodel.o: $(UTILS)transitmodel.f
	$(F77) $(FFLAGS) $(UTILS)transitmodel.f
mandelagol.o: $(UTILS)mandelagol.f
	$(F77) $(FFLAGS) $(UTILS)mandelagol.f
dataio.o: $(UTILS)dataio.f
	$(F77) $(FFLAGS) $(UTILS)dataio.f
fittransitmodel.o: $(UTILS)fittransitmodel.f
	$(F77) $(FFLAGS) $(UTILS)fittransitmodel.f
mrqmin.o: $(UTILS)mrqmin.f
	$(F77) $(FFLAGS) $(UTILS)mrqmin.f
keplerian.o: $(UTILS)keplerian.f
	$(F77) $(FFLAGS) $(UTILS)keplerian.f
mcmc.o: $(UTILS)mcmc.f
	$(F77) $(FFLAGS) $(UTILS)mcmc.f
ran2.o: $(UTILS)ran2.f
	$(F77) $(FFLAGS) $(UTILS)ran2.f
	
	
# Removing object files
.PHONY : clean
clean :
	rm *.o
