#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F77 = ifort
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
# libraries for linking PGPLOT
LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng
#LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng 
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/bin/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: transitfit5 transitmcmc5 mcmchist5 mcmc_corr mcmcphys mcmchistnew a2n plottimeseries rvplot transitplot5 updaten transitchisq5 transitdur5 n1update transitcut5 transitsigclip5 transittiming5 transitsn5

sn5incl = dataio.o getfitpars.o ttcor.o transitmodel.o mandelagol.o keplerian.o occultquad.o stdev.o
transitsn5: transitsn5.f $(sn5incl)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitsn5.f $(sn5incl)

timingincl = dataio.o getfitpars.o ttcor.o transitmodel.o mandelagol.o keplerian.o occultquad.o transitdur.o fittransitmodel.o mrqmin.o
transittiming5: transittiming5.f $(timingincl)
	$(F77) $(LFLAGS) -o $(BIN)$@ transittiming5.f $(timingincl)

prdumpincl = dataio.o getfitpars.o ttcor.o transitmodel.o mandelagol.o keplerian.o occultquad.o
prdump5: prdump5.f $(prdumpincl)
	$(F77) $(LFLAGS) -o $(BIN)$@ prdump5.f $(prdumpincl)

mcmcsetincl = getfitpars.o transitmodel.o mandelagol.o dataio.o fittransitmodel.o mrqmin.o keplerian.o occultquad.o fittransitmodel2.o minpack.o ran2.o getrhosig.o ttcor.o
mcmcsetup: $(mcmcsetincl)
	$(F77) $(LFLAGS) -o $(BIN)$@ mcmcsetup.f $(mcmcsetincl)
tdepth5 = getfitpars.o transitmodel.o mandelagol.o dataio.o fittransitmodel.o mrqmin.o keplerian.o occultquad.o fittransitmodel2.o minpack.o ran2.o getrhosig.o ttcor.o
transitdepth5: $(tdepth5)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitdepth5.f $(tdepth5)
cut5incl = getfitpars.o ttcor.o
transitcut5: $(cut5incl)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitcut5.f $(cut5incl)
occincl = dataio.o getfitpars.o ttcor.o marktransit5.o rqsort.o
occultfind: $(occincl)
	$(F77) $(LFLAGS) -o $(BIN)$@ occultfind.f $(occincl)
n1update: n1update.f getfitpars.o dataio.o
	$(F77) $(LFLAGS) -o $(BIN)$@ n1update.f getfitpars.o dataio.o
datadumpincl = dataio.o getfitpars.o ttcor.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o
datadump5: datadump5.f $(datadumpincl)
	$(F77) $(LFLAGS) -o $(BIN)$@ datadump5.f $(datadumpincl)
durincl = getfitpars.o
transitdur5: transitdur5.f $(durincl)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitdur5.f $(durincl)
remove5incl = getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o ttcor.o
transitremove5: transitremove5.f $(remove5incl)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitremove5.f $(remove5incl)
CHISQINCL = getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o ttcor.o
transitchisq5: transitchisq5.f $(CHISQINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitchisq5.f $(CHISQINCL)
CLIPINCL = getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o ttcor.o
transitsigclip5: transitsigclip5.f $(CLIPINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitsigclip5.f $(CLIPINCL)
SCANINCL5 = getfitpars.o transitmodel.o mandelagol.o dataio.o fittransitmodel.o mrqmin.o keplerian.o occultquad.o fittransitmodel2.o minpack.o ran2.o getrhosig.o ttcor.o
transitscan5: transitscan5.f $(SCANINCL5)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitscan5.f $(SCANINCL5) $(LIBS)
rvplot: rvplot.f
	$(F77) $(LFLAGS) -o $(BIN)$@ rvplot.f $(LIBS)
PTSINCL = dataio.o getfitpars.o
plottimeseries: plottimeseries.f $(PTSINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ plottimeseries.f $(LIBS) $(PTSINCL)
STAR5INCL = dataio.o getfitpars.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o starplot.o
starplot5: starplot5.f $(STAR5INCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ starplot5.f $(LIBS) $(STAR5INCL)
starplot5b: starplot5b.f $(STAR5INCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ starplot5b.f $(LIBS) $(STAR5INCL)
TP5INCL = dataio.o getfitpars.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o
transitplot5: transitplot5.f $(TP5INCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitplot5.f $(LIBS) $(TP5INCL)
UPDATENINCL = getfitpars.o dataio.o
updaten: updaten.f
	$(F77) $(LFLAGS) -o $(BIN)$@ updaten.f $(UPDATENINCL)
a2n: a2n.f
	$(F77) $(LFLAGS) -o $(BIN)$@ a2n.f
HISTNEWINCL = getrhostar.o histogram.o avevar.o rqsort.o getfitpars.o ran2.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o
mcmchistnew: mcmchistnew.f $(HISTNEWINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ mcmchistnew.f $(LIBS) $(HISTNEWINCL)
mcmcphys: mcmcphys.f
	$(F77) $(LFLAGS) -o $(BIN)$@ mcmcphys.f
mcmc_corr: mcmc_corr.f
	$(F77) $(LFLAGS) -o $(BIN)$@ mcmc_corr.f $(LIBS)
HIST5 = getrhostar.o histogram.o avevar.o rqsort.o getfitpars.o ran2.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o
mcmchist5: mcmchist5.f  $(HIST5)
	$(F77) $(LFLAGS) -o $(BIN)$@ mcmchist5.f $(LIBS) $(HIST5)
HIST5bf = getrhostar.o histogram2.o avevar.o rqsort.o getfitpars.o ran2.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o
mcmchist5bf: mcmchist5bf.f  $(HIST5bf)
	$(F77) $(LFLAGS) -o $(BIN)$@ mcmchist5bf.f $(LIBS) $(HIST5bf)
MCMC5INCL = getfitpars.o transitmodel.o mandelagol.o dataio.o fittransitmodel.o mrqmin.o keplerian.o dmcmc.o ran2.o getrhosig.o occultquad.o ttcor.o
transitmcmc5: transitmcmc5.f $(MCMC5INCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitmcmc5.f $(MCMC5INCL)
TRANSINCL3 = getfitpars.o transitmodel.o mandelagol.o dataio.o fittransitmodel.o mrqmin.o keplerian.o occultquad.o fittransitmodel2.o minpack.o getrhosig.o ttcor.o
transitfit5: transitfit5.f $(TRANSINCL3)
	$(F77) $(LFLAGS) -o $(BIN)$@ \
	transitfit5.f $(TRANSINCL3) $(LIBS)
transittiming: transittiming.f $(TRANSINCL3)
	$(F77) $(LFLAGS) -o $(BIN)$@ transittiming.f $(TRANSINCL3) 


#building object libraries
stdev.o: $(UTILS)stdev.f
	$(F77) $(FFLAGS) $(UTILS)stdev.f
transitdur.o: $(UTILS)transitdur.f
	$(F77) $(FFLAGS) $(UTILS)transitdur.f
histogram2.o: $(UTILS)histogram2.f
	$(F77) $(FFLAGS) $(UTILS)histogram2.f
starplot.o: $(UTILS)starplot.f
	$(F77) $(FFLAGS) $(UTILS)starplot.f
marktransit5.o: $(UTILS)marktransit5.f
	$(F77) $(FFLAGS) $(UTILS)marktransit5.f
ttcor.o: $(UTILS)ttcor.f
	$(F77) $(FFLAGS) $(UTILS)ttcor.f
minpack.o: $(UTILS)minpack.f
	$(F77) $(FFLAGS) $(UTILS)minpack.f
occultquad.o: $(UTILS)occultquad.f
	$(F77) $(FFLAGS) $(UTILS)occultquad.f
getrhosig.o: $(UTILS)getrhosig.f
	$(F77) $(FFLAGS) $(UTILS)getrhosig.f
getrhostar.o: $(UTILS)getrhostar.f
	$(F77) $(FFLAGS) $(UTILS)getrhostar.f
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
fittransitmodel2.o: $(UTILS)fittransitmodel2.f
	$(F77) $(FFLAGS) $(UTILS)fittransitmodel2.f
mrqmin.o: $(UTILS)mrqmin.f
	$(F77) $(FFLAGS) $(UTILS)mrqmin.f
keplerian.o: $(UTILS)keplerian.f
	$(F77) $(FFLAGS) $(UTILS)keplerian.f
mcmc.o: $(UTILS)mcmc.f
	$(F77) $(FFLAGS) $(UTILS)mcmc.f
ran2.o: $(UTILS)ran2.f
	$(F77) $(FFLAGS) $(UTILS)ran2.f
dmcmc.o: $(UTILS)dmcmc.f
	$(F77) $(FFLAGS) $(UTILS)dmcmc.f
histogram.o: $(UTILS)histogram.f
	$(F77) $(FFLAGS) $(UTILS)histogram.f
avevar.o: $(UTILS)avevar.f
	$(F77) $(FFLAGS) $(UTILS)avevar.f	
rqsort.o: $(UTILS)rqsort.f
	$(F77) $(FFLAGS) $(UTILS)rqsort.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
