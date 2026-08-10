# ==============================================================================
# Makefile for datatest
# ==============================================================================

FC = gfortran
CC = gcc
OPT ?= -O3
OPENMP ?= -fopenmp
FFLAGS ?= $(OPT) $(OPENMP) -std=legacy -fallow-argument-mismatch
CFLAGS ?= $(OPT) -std=gnu89 -Wno-implicit-function-declaration

# Output directory for binaries
BIN ?= ./bin/

# Utilities subdirectory
UTILS = utils/

# Operating system detection
UNAME_S := $(shell uname -s)

# PGPLOT library detection
ifeq ($(PGPLOT_DIR),)
  ifneq ($(wildcard $(HOME)/Software/pgplot/build/libpgplot.a),)
    PGPLOT_DIR = $(HOME)/Software/pgplot/build
  else ifneq ($(wildcard /Users/jasonrowe/Software/pgplot/build/libpgplot.a),)
    PGPLOT_DIR = /Users/jasonrowe/Software/pgplot/build
  else ifneq ($(wildcard /usr/local/lib/libpgplot.so /usr/local/lib/libpgplot.a),)
    PGPLOT_DIR = /usr/local/lib
  endif
endif

ifneq ($(PGPLOT_DIR),)
  PGPLOT_LFLAGS ?= -L$(PGPLOT_DIR) -lpgplot
endif

ifeq ($(UNAME_S),Darwin)
  X11_LFLAGS ?= -L/opt/homebrew/lib -L/opt/X11/lib -lX11 -lpng -lz
else
  X11_LFLAGS ?= -L/usr/X11R6/lib -L/usr/lib64 -L/usr/local/lib -lX11 -lpng -lz
endif

LIBS_PGPLOT = $(PGPLOT_LFLAGS) $(X11_LFLAGS)

# CFITSIO library detection
CFITSIO_CFLAGS := $(shell pkg-config --cflags cfitsio 2>/dev/null)
CFITSIO_LFLAGS := $(shell pkg-config --libs cfitsio 2>/dev/null)
ifeq ($(CFITSIO_LFLAGS),)
  ifneq ($(wildcard /opt/homebrew/include/fitsio.h),)
    CFITSIO_CFLAGS = -I/opt/homebrew/include
    CFITSIO_LFLAGS = -L/opt/homebrew/lib -lcfitsio
  else ifneq ($(wildcard /usr/local/include/fitsio.h),)
    CFITSIO_CFLAGS = -I/usr/local/include
    CFITSIO_LFLAGS = -L/usr/local/lib -lcfitsio
  else
    CFITSIO_LFLAGS = -lcfitsio
  endif
endif

LIBS_CFITSIO = $(CFITSIO_CFLAGS) $(CFITSIO_LFLAGS)

# Targets to build
TARGETS = detrend51 detrend5 detrend5ramp combine sigclip kfitsread detrend2 \
          bjdcorr dvmatch gmassrad klccomb lcupdateformat ftest transitcut \
          readkeplersc readkeplerlc freqremove plottest kicpick kicstats \
          q1read detrend

.PHONY: all clean init_bin

all: init_bin $(TARGETS)

init_bin:
	@mkdir -p $(BIN)

# Individual binary rules
D51_OBJS = precision.o readdata2.o getfitpars.o ttcor.o marktransit51.o \
           polyfilter51.o rqsort.o lfit.o gaussj.o stdev.o
detrend51: detrend51.f90 $(D51_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ detrend51.f90 $(D51_OBJS)

D5_OBJS = stdev.o rqsort.o polyfilter.o boxcar.o ttcor.o marktransit5.o mrqmin.o
detrend5: detrend5.f $(D5_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ detrend5.f $(D5_OBJS)

D5RAMP_OBJS = stdev.o rqsort.o polyfilter_ramp.o boxcar.o ttcor.o marktransit5.o mrqmin.o
detrend5ramp: detrend5.f $(D5RAMP_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ detrend5.f $(D5RAMP_OBJS)

combine: combine.f90 precision.o rqsort.o
	$(FC) $(FFLAGS) -o $(BIN)$@ combine.f90 precision.o rqsort.o

sigclip: sigclip.f stdev.o marktransit.o transitdur.o
	$(FC) $(FFLAGS) -o $(BIN)$@ sigclip.f stdev.o marktransit.o transitdur.o

kfitsread: kfitsread.f
	$(FC) $(FFLAGS) $(LIBS_CFITSIO) -o $(BIN)$@ kfitsread.f

D2_OBJS = stdev.o rqsort.o polyfilter.o marktransit.o transitdur.o mrqmin.o boxcar.o
detrend2: detrend2.f $(D2_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ detrend2.f $(D2_OBJS)

bjdcorr: bjdcorr.f
	$(FC) $(FFLAGS) -o $(BIN)$@ bjdcorr.f

dvmatch: dvmatch.f
	$(FC) $(FFLAGS) -o $(BIN)$@ dvmatch.f

gmassrad: gmassrad.f
	$(FC) $(FFLAGS) -o $(BIN)$@ gmassrad.f

klccomb: klccomb.f
	$(FC) $(FFLAGS) -o $(BIN)$@ klccomb.f

lcupdateformat: lcupdateformat.f90 precision.o
	$(FC) $(FFLAGS) -o $(BIN)$@ lcupdateformat.f90 precision.o

ftest: ftest.f
	$(FC) $(FFLAGS) -o $(BIN)$@ ftest.f

transitcut: transitcut.f marktransit.o transitdur.o
	$(FC) $(FFLAGS) -o $(BIN)$@ transitcut.f marktransit.o transitdur.o

readkeplersc: readkeplersc.f
	$(FC) $(FFLAGS) -o $(BIN)$@ readkeplersc.f $(LIBS_PGPLOT)

readkeplerlc: readkeplerlc.f
	$(FC) $(FFLAGS) -o $(BIN)$@ readkeplerlc.f $(LIBS_PGPLOT)

freqremove: freqremove.f
	$(FC) $(FFLAGS) -o $(BIN)$@ freqremove.f $(LIBS_PGPLOT)

plottest: plottest.f
	$(FC) $(FFLAGS) -o $(BIN)$@ plottest.f $(LIBS_PGPLOT)

kicpick: kicpick.f
	$(FC) $(FFLAGS) -o $(BIN)$@ kicpick.f

kicstats: kicstats.f
	$(FC) $(FFLAGS) -o $(BIN)$@ kicstats.f

q1read: q1read.f
	$(FC) $(FFLAGS) -o $(BIN)$@ q1read.f $(LIBS_PGPLOT)

DETREND_OBJS = boxcar.o dataio.o stdev.o polydetrend.o mrqmin.o rqsort.o
detrend: detrend.f $(DETREND_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ detrend.f $(DETREND_OBJS)

# Object compilation rules
precision.o: $(UTILS)precision.f90
	$(FC) -c $(FFLAGS) -o $@ $<

readdata2.o: $(UTILS)readdata2.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

marktransit51.o: $(UTILS)marktransit51.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

polyfilter51.o: $(UTILS)polyfilter51.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

%.o: $(UTILS)%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%.o: $(UTILS)%.f
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod *__genmod.f90 $(BIN)detrend51 $(BIN)detrend5 \
	      $(BIN)detrend5ramp $(BIN)combine $(BIN)sigclip $(BIN)kfitsread \
	      $(BIN)detrend2 $(BIN)bjdcorr $(BIN)dvmatch $(BIN)gmassrad \
	      $(BIN)klccomb $(BIN)lcupdateformat $(BIN)ftest $(BIN)transitcut \
	      $(BIN)readkeplersc $(BIN)readkeplerlc $(BIN)freqremove \
	      $(BIN)plottest $(BIN)kicpick $(BIN)kicstats $(BIN)q1read $(BIN)detrend
