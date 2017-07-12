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

all: transitfit2 modelbuild transitboot transitplot transiterr cadextract transitfit3 transitboot2.0 transitmcmc3 transitdur prdump transitremove

claretquad: claretquad.f
	$(F77) $(LFLAGS) -o $(BIN)$@ claretquad.f

claretlimb: claretlimb.f
	$(F77) $(LFLAGS) -o $(BIN)$@ claretlimb.f

PLOT3INCL = dataio.o getfitpars.o transitmodel.o keplerian.o mandelagol.o eclmod.o mathematics.o getplotpars.o mag2flux.o graphics.o plotflux.o plotphflux.o modelplotflux.o binp.o plotoddeven.o transn.o
transitplot3: transitplot3.f $(PLOT3INCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitplot3.f $(PLOT3INCL) $(LIBS)

PLOT2INCL = dataio.o getfitpars.o transitmodel.o keplerian.o mandelagol.o eclmod.o mathematics.o
transitplot2: transitplot2.f $(PLOT2INCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitplot2.f $(PLOT2INCL) $(LIBS)

TREMOVE = getfitpars.o dataio.o transitmodel.o keplerian.o mandelagol.o eclmod.o
transitremove: transitremove.f $(TREMOVE)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitremove.f $(TREMOVE)

transitadd: transitadd.f $(TREMOVE)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitadd.f $(TREMOVE)

TSNINCL = getfitpars.o transitmodel.o keplerian.o mandelagol.o eclmod.o dataio.o
transitsn: transitsn.f	$(TSNINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitsn.f $(TSNINCL)

TDEPINCL = getfitpars.o transitmodel.o keplerian.o mandelagol.o eclmod.o
transitdepth: transitdepth.f $(TDEPINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitdepth.f $(TDEPINCL)

DENINCL = getfitpars.o
densityout: densityout.f $(DENINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ densityout.f $(DENINCL)

DUMPINCL = dataio.o getfitpars.o keplerian.o mandelagol.o eclmod.o
datadump: datadump.f $(DUMPINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ datadump.f $(DUMPINCL)

PRINCL = dataio.o getfitpars.o keplerian.o mandelagol.o eclmod.o
prdump: prdump.f $(PRINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ prdump.f $(PRINCL)

LIMBINCL = mandelagolv2.o gamma.o f21r.o isneg.o get_transformation.o check_parms.o f1conv.o f2.o f2gam.o hypgeof1.o writef1.o f1bnl.o hypgeo.o hypstartf1.o
limbtest: limbtest.f $(LIMBINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ limbtest.f $(LIMBINCL)

transitdur: transitdur.f getfitpars.o
	$(F77) $(LFLAGS) -o $(BIN)$@ transitdur.f getfitpars.o

updatef: updatef.f getfitpars.o dataio.o
	$(F77) $(LFLAGS) -o $(BIN)$@ updatef.f getfitpars.o dataio.o

MCMCINCL3 = dataio.o getfitpars.o mathematics.o transitmodel.o keplerian.o mandelagol.o eclmod.o mcmc.o
transitmcmc3: transitmcmc3.f $(MCMCINCL3)
	$(F77) $(LFLAGS) -o $(BIN)$@ transitmcmc3.f $(MCMCINCL3)

boothist: boothist.f
	$(F77) $(LFLAGS) -o $(BIN)$@ boothist.f $(LIBS)

adrs: adrs.f
	$(F77) $(LFLAGS) -o $(BIN)$@ adrs.f

limbdarken: limbdarken.f
	$(F77) $(LFLAGS) -o $(BIN)$@ limbdarken.f

CADINCL = dataio.o getfitpars.o
cadextract: cadextract.f $(CADINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ cadextract.f $(CADINCL)

tps_sort: tps_sort.f
	$(F77) $(LFLAGS) -o $(BIN)$@ tps_sort.f

TPSINCL = dataio.o boxcar.o stdev.o rqsort.o
tpsguess: tpsguess.f $(TPSINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ tpsguess.f $(TPSINCL)

tpsguess2: tpsguess2.f $(TPSINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ tpsguess2.f $(TPSINCL)

TESTINCL = dataio.o transitmodel.o keplerian.o mandelagol.o eclmod.o
testfit: testfit.f $(TESTINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ testfit.f $(ERRINCL)

ERRINCL = dataio.o getfitpars.o marginalize.o likelihood.o transitmodel.o keplerian.o mandelagol.o eclmod.o
transiterr: transiterr.f $(ERRINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ \
	transiterr.f $(ERRINCL) 

PLOTINCL = dataio.o getfitpars.o transitmodel.o keplerian.o mandelagol.o eclmod.o mathematics.o getplotpars.o mag2flux.o binp.o
transitplot: transitplot.f $(PLOTINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ \
	transitplot.f $(PLOTINCL) $(LIBS)

TRANSINCL = getfitpars.o transitmodel.o modelplot.o eclmod.o mandelagol.o dataio.o graphics.o mathematics.o fittransitmodel.o cfunk.o likelihood.o amoeba.o autoexplore.o keplerian.o
transitfit2: transitfit2.0.f $(TRANSINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ \
	transitfit2.0.f $(TRANSINCL) $(LIBS)
TRANSINCL3 = getfitpars.o transitmodel.o modelplot.o eclmod.o mandelagol.o dataio.o graphics.o mathematics.o fittransitmodel2.o likelihood.o autoexplore.o keplerian.o mrqmin.o dfridr.o boxcar.o oddeven.o
transitfit3: transitfit3.0.f $(TRANSINCL3)
	$(F77) $(LFLAGS) -o $(BIN)$@ \
	transitfit3.0.f $(TRANSINCL3) $(LIBS)

BUILDINCL = dataio.o graphics.o getfitpars.o mathematics.o modelplot.o transitmodel.o eclmod.o mandelagol.o explore.o likelihood.o keplerian.o
modelbuild: modelbuild.f $(BUILDINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ \
	modelbuild.f $(BUILDINCL) $(LIBS)

BOOTINCL = getfitpars.o transitmodel.o eclmod.o mandelagol.o dataio.o mathematics.o fittransitmodel2.o likelihood.o keplerian.o mrqmin.o
transitboot: transitboot.f $(BOOTINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ \
	transitboot.f $(BOOTINCL)

transitboot2.0: transitboot2.0.f $(BOOTINCL)
	$(F77) $(LFLAGS) -o $(BIN)$@ \
	transitboot2.0.f $(BOOTINCL)

#building object libraries

rqsort.o: $(UTILS)rqsort.f
	$(F77) $(FFLAGS) $(UTILS)rqsort.f
transn.o: $(UTILS)transn.f
	$(F77) $(FFLAGS) $(UTILS)transn.f
plotoddeven.o: $(UTILS)plotoddeven.f
	$(F77) $(FFLAGS) $(UTILS)plotoddeven.f
binp.o: $(UTILS)binp.f
	$(F77) $(FFLAGS) $(UTILS)binp.f
modelplotflux.o: $(UTILS)modelplotflux.f
	$(F77) $(FFLAGS) $(UTILS)modelplotflux.f
plotphflux.o: $(UTILS)plotphflux.f
	$(F77) $(FFLAGS) $(UTILS)plotphflux.f
plotflux.o: $(UTILS)plotflux.f
	$(F77) $(FFLAGS) $(UTILS)plotflux.f
mag2flux.o: $(UTILS)mag2flux.f	
	$(F77) $(FFLAGS) $(UTILS)mag2flux.f
getplotpars.o: $(UTILS)getplotpars.f
	$(F77) $(FFLAGS) $(UTILS)getplotpars.f
hypstartf1.o: $(UTILS)hypstartf1.f
	$(F77) $(FFLAGS) $(UTILS)hypstartf1.f
hypgeo.o: $(UTILS)hypgeo.f
	$(F77) $(FFLAGS) $(UTILS)hypgeo.f
f1bnl.o: $(UTILS)f1bnl.f
	$(F77) $(FFLAGS) $(UTILS)f1bnl.f
writef1.o: $(UTILS)writef1.f
	$(F77) $(FFLAGS) $(UTILS)writef1.f
hypgeof1.o: $(UTILS)hypgeof1.f
	$(F77) $(FFLAGS) $(UTILS)hypgeof1.f
f2gam.o: $(UTILS)f2gam.f
	$(F77) $(FFLAGS) $(UTILS)f2gam.f
f2.o: $(UTILS)f2.f
	$(F77) $(FFLAGS) $(UTILS)f2.f
f1conv.o: $(UTILS)f1conv.f
	$(F77) $(FFLAGS) $(UTILS)f1conv.f
check_parms.o: $(UTILS)check_parms.f
	$(F77) $(FFLAGS) $(UTILS)check_parms.f
get_transformation.o: $(UTILS)get_transformation.f
	$(F77) $(FFLAGS) $(UTILS)get_transformation.f
isneg.o: $(UTILS)isneg.f
	$(F77) $(FFLAGS) $(UTILS)isneg.f
f21r.o: $(UTILS)f21r.f
	$(F77) $(FFLAGS) $(UTILS)f21r.f
gamma.o: $(UTILS)gamma.f
	$(F77) $(FFLAGS) $(UTILS)gamma.f
mandelagolv2.o: $(UTILS)mandelagolv2.f
	$(F77) $(FFLAGS) $(UTILS)mandelagolv2.f
stdev.o: $(UTILS)stdev.f
	$(F77) $(FFLAGS) $(UTILS)stdev.f
getfitpars.o: $(UTILS)getfitpars.f
	$(F77) $(FFLAGS) $(UTILS)getfitpars.f
transitmodel.o: $(UTILS)transitmodel.f
	$(F77) $(FFLAGS) $(UTILS)transitmodel.f
modelplot.o: $(UTILS)modelplot.f
	$(F77) $(FFLAGS) $(UTILS)modelplot.f
eclmod.o: $(UTILS)eclmod.f
	$(F77) $(FFLAGS) $(UTILS)eclmod.f
mandelagol.o: $(UTILS)mandelagol.f
	$(F77) $(FFLAGS) $(UTILS)mandelagol.f
dataio.o: $(UTILS)dataio.f
	$(F77) $(FFLAGS) $(UTILS)dataio.f
graphics.o: $(UTILS)graphics.f
	$(F77) $(FFLAGS) $(UTILS)graphics.f
mathematics.o: $(UTILS)mathematics.f
	$(F77) $(FFLAGS) $(UTILS)mathematics.f
fittransitmodel.o: $(UTILS)fittransitmodel.f
	$(F77) $(FFLAGS) $(UTILS)fittransitmodel.f
cfunk.o: $(UTILS)cfunk.f
	$(F77) $(FFLAGS) $(UTILS)cfunk.f
amoeba.o: $(UTILS)amoeba.f
	$(F77) $(FFLAGS) $(UTILS)amoeba.f
explore.o: $(UTILS)explore.f
	$(F77) $(FFLAGS) $(UTILS)explore.f
autoexplore.o: $(UTILS)autoexplore.f
	$(F77) $(FFLAGS) $(UTILS)autoexplore.f
keplerian.o: $(UTILS)keplerian.f
	$(F77) $(FFLAGS) $(UTILS)keplerian.f
marginalize.o: $(UTILS)marginalize.f
	$(F77) $(FFLAGS) $(UTILS)marginalize.f
likelihood.o: $(UTILS)likelihood.f
	$(F77) $(FFLAGS) $(UTILS)likelihood.f
fittransitmodel2.o: $(UTILS)fittransitmodel2.f
	$(F77) $(FFLAGS) $(UTILS)fittransitmodel2.f
mrqmin.o: $(UTILS)mrqmin.f
	$(F77) $(FFLAGS) $(UTILS)mrqmin.f
dfridr.o: $(UTILS)dfridr.f
	$(F77) $(FFLAGS) $(UTILS)dfridr.f
boxcar.o: $(UTILS)boxcar.f
	$(F77) $(FFLAGS) $(UTILS)boxcar.f
oddeven.o: $(UTILS)oddeven.f
	$(F77) $(FFLAGS) $(UTILS)oddeven.f
mcmc.o: $(UTILS)mcmc.f
	$(F77) $(FFLAGS) $(UTILS)mcmc.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
