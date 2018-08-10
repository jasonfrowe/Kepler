#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F77 = ifort
F90 = ifort
#compiling object file flags
OPT1 = -O3 
OPT2 = 
FFLAGS = -c $(OPT1) $(OPT2)
#FFLAGS = -c -O0 -g -CB -warn $(OPT2)
#linking flags
LFLAGS = $(OPT1) $(OPT2)
#LFLAGS = -O0 -g -CB -warn $(OPT2)
#other linking flags
#LFLAGS2 = -O2
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

all: keplerper stdplot ttper freqmod

freqmodincl = readdata.o precision.o nyquest.o rqsort.o sigclip.o stdev.o bindt.o avgrm.o moment.o graphics.o jmfourw.o phasept.o cosinefit.o shaperm.o windowfn.o mrqmin.o integrate.o readfreqs.o freqout.o fixa.o
freqmod: freqmod.f90 $(freqmodincl)
	$(F90) $(LFLAGS) -o /Users/rowe/Documents/Kepler/freqmod/$@ freqmod.f90 $(freqmodincl) $(LIBS)

perincls = readdata.o nyquest.o modecal.o stdev.o sigclip.o avgrm.o moment.o bindt.o graphics.o jmfourw.o ovrwrt.o cosinefit.o integrate.o phasept.o rqsort.o mrqmin.o windowfn.o shaperm.o detrend.o boxcar_old.o plotfit.o freqout.o
keplerper: keplerper.f $(perincls)
	$(F77) $(LFLAGS) -o $(BIN)$@ keplerper.f $(perincls) $(LIBS)

ttincls = readdata_tt.o nyquest.o modecal.o stdev.o sigclip.o avgrm.o moment.o bindt.o graphics_tt.o jmfourw_tt.o ovrwrt.o cosinefit.o integrate.o phasept.o rqsort.o mrqmin.o windowfn.o shaperm.o detrend.o boxcar_old.o plotfit.o freqout.o
ttper: ttper.f $(ttincls)
	$(F77) $(LFLAGS) -o $(BIN)$@ ttper.f $(ttincls) $(LIBS)

quickplot: quickplot.f $(perincls) sigup.o
	$(F77) $(LFLAGS) -o $(BIN)$@ quickplot.f $(perincls) sigup.o $(LIBS)

stdplot: stdplot.f
	$(F77) $(LFLAGS) -o $(BIN)$@ stdplot.f $(LIBS)

#building object libraries

precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
freqout.o: $(UTILS)freqout.f
	$(F77) $(FFLAGS) $(UTILS)freqout.f
plotfit.o: $(UTILS)plotfit.f
	$(F77) $(FFLAGS) $(UTILS)plotfit.f
readdata.o: $(UTILS)readdata.f
	$(F77) $(FFLAGS) $(UTILS)readdata.f
readdata_tt.o: $(UTILS)readdata_tt.f
	$(F77) $(FFLAGS) $(UTILS)readdata_tt.f
nyquest.o: $(UTILS)nyquest.f
	$(F77) $(FFLAGS) $(UTILS)nyquest.f
modecal.o: $(UTILS)modecal.f
	$(F77) $(FFLAGS) $(UTILS)modecal.f
stdev.o: $(UTILS)stdev.f
	$(F77) $(FFLAGS) $(UTILS)stdev.f
sigclip.o: $(UTILS)sigclip.f
	$(F77) $(FFLAGS) $(UTILS)sigclip.f
avgrm.o: $(UTILS)avgrm.f
	$(F77) $(FFLAGS) $(UTILS)avgrm.f
moment.o: $(UTILS)moment.f
	$(F77) $(FFLAGS) $(UTILS)moment.f
bindt.o: $(UTILS)bindt.f
	$(F77) $(FFLAGS) $(UTILS)bindt.f
graphics.o: $(UTILS)graphics.f
	$(F77) $(FFLAGS) $(UTILS)graphics.f
graphics_tt.o: $(UTILS)graphics_tt.f
	$(F77) $(FFLAGS) $(UTILS)graphics_tt.f
jmfourw.o: $(UTILS)jmfourw.f
	$(F77) $(FFLAGS) $(UTILS)jmfourw.f
jmfourw_tt.o: $(UTILS)jmfourw_tt.f
	$(F77) $(FFLAGS) $(UTILS)jmfourw_tt.f
ovrwrt.o: $(UTILS)ovrwrt.f
	$(F77) $(FFLAGS) $(UTILS)ovrwrt.f
cosinefit.o: $(UTILS)cosinefit.f
	$(F77) $(FFLAGS) $(UTILS)cosinefit.f
integrate.o: $(UTILS)integrate.f
	$(F77) $(FFLAGS) $(UTILS)integrate.f
phasept.o: $(UTILS)phasept.f
	$(F77) $(FFLAGS) $(UTILS)phasept.f
rqsort.o: $(UTILS)rqsort.f
	$(F77) $(FFLAGS) $(UTILS)rqsort.f
mrqmin.o: $(UTILS)mrqmin.f
	$(F77) $(FFLAGS) $(UTILS)mrqmin.f
windowfn.o: $(UTILS)windowfn.f
	$(F77) $(FFLAGS) $(UTILS)windowfn.f
shaperm.o: $(UTILS)shaperm.f
	$(F77) $(FFLAGS) $(UTILS)shaperm.f
detrend.o: $(UTILS)detrend.f
	$(F77) $(FFLAGS) $(UTILS)detrend.f
boxcar.o: $(UTILS)boxcar.f
	$(F77) $(FFLAGS) $(UTILS)boxcar.f
boxcar_old.o: $(UTILS)boxcar_old.f
	$(F77) $(FFLAGS) $(UTILS)boxcar_old.f
sigup.o: $(UTILS)sigup.f
	$(F77) $(FFLAGS) $(UTILS)sigup.f
readfreqs.o: $(UTILS)readfreqs.f
	$(F77) $(FFLAGS) $(UTILS)readfreqs.f
fixa.o: $(UTILS)fixa.f
	$(F77) $(FFLAGS) $(UTILS)fixa.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
