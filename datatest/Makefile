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
#other linking flags
LFLAGS2 = -O2
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

all: detrend2 kfitsread detrend5 sigclip freqremove bjdcorr

plottest: plottest.f
	$(F77) $(LFLAGS) -o /home/rowe/temp/$@ plottest.f $(LIBS)
kicpick: kicpick.f
	$(F77) $(LFLAGS) -o /home/rowe/Kepler/kicwork/$@ kicpick.f
kicstats: kicstats.f
	$(F77) $(LFLAGS) -o /home/rowe/Kepler/kicwork/$@ kicstats.f
kfitsread: kfitsread.f
	$(F77) $(LFLAGS) -o$(BIN)$@ kfitsread.f $(LIBS2)
detrend5incl = stdev.o rqsort.o polyfilter.o boxcar.o ttcor.o marktransit5.o mrqmin.o
detrend5: detrend5.f $(detrend5incl)
	$(F77) $(LFLAGS) -o$(BIN)$@ detrend5.f $(detrend5incl)
ftest: ftest.f
	$(F77) $(LFLAGS) -o$(BIN)$@ ftest.f
transitcutincl = marktransit.o transitdur.o
transitcut: transitcut.f $(transitcutincl)
	$(F77) $(LFLAGS) -o$(BIN)$@ transitcut.f $(transitcutincl)
sigclipincl = stdev.o marktransit.o transitdur.o
sigclip: sigclip.f $(sigclipincl)
	$(F77) $(LFLAGS) -o$(BIN)$@ sigclip.f $(sigclipincl)
detrend2incl = stdev.o rqsort.o polyfilter.o marktransit.o transitdur.o mrqmin.o boxcar.o
detrend2: detrend2.f $(detrend2incl)
	$(F77) $(LFLAGS) -o $(BIN)$@ detrend2.f $(detrend2incl)
dvmatch: dvmatch.f
	$(F77) $(LFLAGS) -o $(BIN)$@ dvmatch.f
bjdcorr: bjdcorr.f
	$(F77) $(LFLAGS) -o $(BIN)$@ bjdcorr.f
gmassrad: gmassrad.f
	$(F77) $(LFLAGS) -o $(BIN)$@ gmassrad.f
klccomb: klccomb.f
	$(F77) $(LFLAGS) -o $(BIN)$@ klccomb.f
q1read: q1read.f
	$(F77) $(LFLAGS) -o $(BIN)$@ q1read.f $(LIBS)
freqremove: freqremove.f
	$(F77) $(LFLAGS) -o $(BIN)$@ freqremove.f $(LIBS)
readkeplersc: readkeplersc.f
	$(F77) $(LFLAGS) -o $(BIN)$@ readkeplersc.f $(LIBS)
readkeplerlc: readkeplerlc.f 
	$(F77) $(LFLAGS) -o $(BIN)$@ readkeplerlc.f $(LIBS)
detrendincl = boxcar.o dataio.o stdev.o polydetrend.o mrqmin.o rqsort.o
detrend: detrend.f $(detrendincl)
	$(F77) $(LFLAGS) -o $(BIN)$@ detrend.f $(detrendincl)

#building object libraries

marktransit5.o: $(UTILS)marktransit5.f
	$(F77) $(FFLAGS) $(UTILS)marktransit5.f
ttcor.o: $(UTILS)ttcor.f
	$(F77) $(FFLAGS) $(UTILS)ttcor.f
transitdur.o: $(UTILS)transitdur.f
	$(F77) $(FFLAGS) $(UTILS)transitdur.f
marktransit.o: $(UTILS)marktransit.f
	$(F77) $(FFLAGS) $(UTILS)marktransit.f
polyfilter.o: $(UTILS)polyfilter.f
	$(F77) $(FFLAGS) $(UTILS)polyfilter.f
rqsort.o: $(UTILS)rqsort.f
	$(F77) $(FFLAGS) $(UTILS)rqsort.f
mrqmin.o: $(UTILS)mrqmin.f
	$(F77) $(FFLAGS) $(UTILS)mrqmin.f
polydetrend.o: $(UTILS)polydetrend.f
	$(F77) $(FFLAGS) $(UTILS)polydetrend.f
boxcar.o: $(UTILS)boxcar.f
	$(F77) $(FFLAGS) $(UTILS)boxcar.f
dataio.o: $(UTILS)dataio.f
	$(F77) $(FFLAGS) $(UTILS)dataio.f
stdev.o: $(UTILS)stdev.f
	$(F77) $(FFLAGS) $(UTILS)stdev.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
