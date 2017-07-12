#Name of C compiler
CCMPL = icc
#Name of Fortran compiler
F77 = ifort
#F77 = gfortran
#compiling object file flags
OPT1 = -O3 -ftz -ip -ipo -march=core2 -mtune=core2 -assume buffered_io -warn
OPT2 = -warn alignments -align all -align rec8byte -warn -mcmodel medium -shared-intel
#OPT1 = -O3 -mcmodel medium -shared-intel
#OPT2 = 
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
#LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng 
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng -lsvml
#Directory where executable are placed
BIN = /home/rowe/p555/transitfit/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: transitfit6 transitmcmc6 transitplot6 mcmchist6

HISTINCL6 = getfitpars.o ran2.o histogram.o getrhostar.o transitmodel.o keplerian.o avevar.o rqsort.o
mcmchist6: mcmchist6.f $(HISTINCL6)
	$(F77) $(LFLAGS) -o $(BIN)$@ mcmchist6.f $(HISTINCL6) $(LIBS)
PLOTINCL6 = getfitpars.o dataio.o keplerian.o transitdur.o graphics.o transitmodel.o
transitplot6: transitplot6.f $(PLOTINCL6)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitplot6.f $(PLOTINCL6) $(LIBS)	
TRANSINCL6 = getfitpars.o dataio.o keplerian.o transitdur.o fittransitmodel2.o minpack.o graphics.o transitmodel.o
transitfit6: transitfit6.f $(TRANSINCL6)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitfit6.f $(TRANSINCL6) $(LIBS)
MCMCINCL6 = getfitpars.o dataio.o keplerian.o transitdur.o graphics.o transitmodel.o ran2.o dmcmc.o getrhosig.o
transitmcmc6: transitmcmc6.f $(MCMCINCL6)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitmcmc6.f $(MCMCINCL6) $(LIBS)

#building object libraries
avevar.o: $(UTILS)avevar.f
	$(F77) $(FFLAGS) $(UTILS)avevar.f
rqsort.o: $(UTILS)rqsort.f
	$(F77) $(FFLAGS) $(UTILS)rqsort.f
getrhostar.o: $(UTILS)getrhostar.f
	$(F77) $(FFLAGS) $(UTILS)getrhostar.f
histogram.o: $(UTILS)histogram.f
	$(F77) $(FFLAGS) $(UTILS)histogram.f
getrhosig.o: $(UTILS)getrhosig.f
	$(F77) $(FFLAGS) $(UTILS)getrhosig.f
dmcmc.o: $(UTILS)dmcmc.f
	$(F77) $(FFLAGS) $(UTILS)dmcmc.f
ran2.o: $(UTILS)ran2.f
	$(F77) $(FFLAGS) $(UTILS)ran2.f
transitmodel.o: $(UTILS)transitmodel.f
	$(F77) $(FFLAGS) $(UTILS)transitmodel.f
graphics.o: $(UTILS)graphics.f
	$(F77) $(FFLAGS) $(UTILS)graphics.f
minpack.o: $(UTILS)minpack.f
	$(F77) $(FFLAGS) $(UTILS)minpack.f
mrqmin.o: $(UTILS)mrqmin.f
	$(F77) $(FFLAGS) $(UTILS)mrqmin.f
fittransitmodel2.o: $(UTILS)fittransitmodel2.f
	$(F77) $(FFLAGS) $(UTILS)fittransitmodel2.f
fittransitmodel.o: $(UTILS)fittransitmodel.f
	$(F77) $(FFLAGS) $(UTILS)fittransitmodel.f
keplerian.o: $(UTILS)keplerian.f
	$(F77) $(FFLAGS) $(UTILS)keplerian.f
dataio.o: $(UTILS)dataio.f
	$(F77) $(FFLAGS) $(UTILS)dataio.f
getfitpars.o: $(UTILS)getfitpars.f
	$(F77) $(FFLAGS) $(UTILS)getfitpars.f
transitdur.o: $(UTILS)transitdur.f
	$(F77) $(FFLAGS) $(UTILS)transitdur.f
	
# Removing object files
.PHONY : clean
clean :
	rm *.o
