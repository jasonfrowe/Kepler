# ==============================================================================
# Makefile for transitfit5
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

# Compiler definition macro for CPP preprocessing
ifeq ($(findstring gfortran,$(FC)),gfortran)
  FC_DEF = gfortran
else
  FC_DEF = ifort
endif

# Targets to build
TARGETS = transitfit51 transitfit5 transitmcmc5 transitcut5 sigclip claretquad \
          limbpriors dataclip transitsn5 transittiming5 prdump5 mcmcsetup \
          transitdepth5 occultfind n1update datadump5 transitdur5 transitremove5 \
          transitadd5 transitchisq5 transitsigclip5 transitscan5 rvplot \
          plottimeseries starplot5 starplot5b transitplot5 updaten a2n \
          mcmchistnew mcmcphys mcmc_corr mcmchist5 mcmchist5bf

.PHONY: all clean init_bin

all: init_bin $(TARGETS)

init_bin:
	@mkdir -p $(BIN)

# Common object dependencies
UTILS_OBJS = precision.o fittermod.o getfitpars.o ttcor.o transitmodel.o \
             mandelagol.o keplerian.o occultquad.o fittransitmodel.o \
             fittransitmodel2.o fittransitmodel3.o minpack.o getrhosig.o \
             occultnlbm.o dataio.o dataioSC.o mrqmin.o dmcmc.o ran2.o \
             stdev.o transitdur.o readdata2.o cutoutliers.o meddiff.o \
             calcldprior.o locate.o trilinear.o getrhostar.o histogram.o \
             histogram2.o avevar.o rqsort.o starplot.o marktransit5.o

# Individual binary rules
TF51_OBJS = precision.o fittermod.o getfitpars.o ttcor.o transitmodel.o \
            mandelagol.o keplerian.o occultquad.o fittransitmodel3.o minpack.o \
            getrhosig.o occultnlbm.o
transitfit51: transitfit51.f90 $(TF51_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ transitfit51.f90 $(TF51_OBJS)

TF5_OBJS = getfitpars.o transitmodel.o mandelagol.o dataio.o fittransitmodel.o \
           mrqmin.o keplerian.o occultquad.o fittransitmodel2.o minpack.o \
           getrhosig.o ttcor.o
transitfit5: transitfit5.f $(TF5_OBJS)
	$(FC) -cpp -D$(FC_DEF) $(FFLAGS) -o $(BIN)$@ transitfit5.f $(TF5_OBJS)

MCMC5_OBJS = getfitpars.o transitmodel.o mandelagol.o dataio.o fittransitmodel.o \
             mrqmin.o keplerian.o dmcmc.o ran2.o getrhosig.o occultquad.o ttcor.o
transitmcmc5: transitmcmc5.f $(MCMC5_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ transitmcmc5.f $(MCMC5_OBJS)

transitcut5: transitcut5.f getfitpars.o ttcor.o dataio.o
	$(FC) $(FFLAGS) -o $(BIN)$@ transitcut5.f getfitpars.o ttcor.o dataio.o

SIGCLIP_OBJS = precision.o readdata2.o getfitpars.o ttcor.o stdev.o \
               transitmodel.o mandelagol.o keplerian.o occultquad.o
sigclip: sigclip.f90 $(SIGCLIP_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ sigclip.f90 $(SIGCLIP_OBJS)

claretquad: claretquad.f
	$(FC) $(FFLAGS) -o $(BIN)$@ claretquad.f

limbpriors: limbpriors.f90 precision.o calcldprior.o locate.o trilinear.o
	$(FC) $(FFLAGS) -o $(BIN)$@ limbpriors.f90 precision.o calcldprior.o locate.o trilinear.o

dataclip: dataclip.f90 precision.o cutoutliers.o meddiff.o rqsort.o
	$(FC) $(FFLAGS) -o $(BIN)$@ dataclip.f90 precision.o cutoutliers.o meddiff.o rqsort.o

SN5_OBJS = dataio.o getfitpars.o ttcor.o transitmodel.o mandelagol.o keplerian.o occultquad.o stdev.o
transitsn5: transitsn5.f $(SN5_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ transitsn5.f $(SN5_OBJS)

TIMING_OBJS = dataio.o getfitpars.o ttcor.o transitmodel.o mandelagol.o keplerian.o occultquad.o transitdur.o fittransitmodel.o mrqmin.o
transittiming5: transittiming5.f $(TIMING_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ transittiming5.f $(TIMING_OBJS)

prdump5: prdump5.f dataio.o getfitpars.o ttcor.o transitmodel.o mandelagol.o keplerian.o occultquad.o
	$(FC) $(FFLAGS) -o $(BIN)$@ prdump5.f dataio.o getfitpars.o ttcor.o transitmodel.o mandelagol.o keplerian.o occultquad.o

MCMCSET_OBJS = getfitpars.o transitmodel.o mandelagol.o dataio.o fittransitmodel.o mrqmin.o keplerian.o occultquad.o fittransitmodel2.o minpack.o ran2.o getrhosig.o ttcor.o
mcmcsetup: mcmcsetup.f $(MCMCSET_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ mcmcsetup.f $(MCMCSET_OBJS)

transitdepth5: transitdepth5.f $(MCMCSET_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ transitdepth5.f $(MCMCSET_OBJS)

occultfind: occultfind.f dataio.o getfitpars.o ttcor.o marktransit5.o rqsort.o
	$(FC) $(FFLAGS) -o $(BIN)$@ occultfind.f dataio.o getfitpars.o ttcor.o marktransit5.o rqsort.o

n1update: n1update.f getfitpars.o dataio.o
	$(FC) $(FFLAGS) -o $(BIN)$@ n1update.f getfitpars.o dataio.o

datadump5: datadump5.f getfitpars.o ttcor.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o
	$(FC) $(FFLAGS) -o $(BIN)$@ datadump5.f getfitpars.o ttcor.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o

transitdur5: transitdur5.f getfitpars.o
	$(FC) $(FFLAGS) -o $(BIN)$@ transitdur5.f getfitpars.o

transitremove5: transitremove5.f getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o ttcor.o
	$(FC) $(FFLAGS) -o $(BIN)$@ transitremove5.f getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o ttcor.o

transitadd5: transitadd5.f getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o ttcor.o
	$(FC) $(FFLAGS) -o $(BIN)$@ transitadd5.f getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o ttcor.o

transitchisq5: transitchisq5.f getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o ttcor.o
	$(FC) $(FFLAGS) -o $(BIN)$@ transitchisq5.f getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o ttcor.o

transitsigclip5: transitsigclip5.f getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o ttcor.o
	$(FC) $(FFLAGS) -o $(BIN)$@ transitsigclip5.f getfitpars.o transitmodel.o mandelagol.o dataio.o keplerian.o occultquad.o ttcor.o

transitscan5: transitscan5.f $(MCMCSET_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ transitscan5.f $(MCMCSET_OBJS) $(LIBS_PGPLOT)

rvplot: rvplot.f
	$(FC) $(FFLAGS) -o $(BIN)$@ rvplot.f $(LIBS_PGPLOT)

plottimeseries: plottimeseries.f dataio.o getfitpars.o
	$(FC) $(FFLAGS) -o $(BIN)$@ plottimeseries.f dataio.o getfitpars.o $(LIBS_PGPLOT)

starplot5: starplot5.f dataio.o getfitpars.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o starplot.o
	$(FC) $(FFLAGS) -o $(BIN)$@ starplot5.f dataio.o getfitpars.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o starplot.o $(LIBS_PGPLOT)

starplot5b: starplot5b.f dataio.o getfitpars.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o starplot.o
	$(FC) $(FFLAGS) -o $(BIN)$@ starplot5b.f dataio.o getfitpars.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o starplot.o $(LIBS_PGPLOT)

TP5_OBJS = dataio.o getfitpars.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o occultnlbm.o
transitplot5: transitplot5.f $(TP5_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ transitplot5.f $(TP5_OBJS) $(LIBS_PGPLOT)

updaten: updaten.f getfitpars.o dataio.o
	$(FC) $(FFLAGS) -o $(BIN)$@ updaten.f getfitpars.o dataio.o

a2n: a2n.f
	$(FC) $(FFLAGS) -o $(BIN)$@ a2n.f

HIST_OBJS = getrhostar.o histogram.o avevar.o rqsort.o getfitpars.o ran2.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o
mcmchistnew: mcmchistnew.f $(HIST_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ mcmchistnew.f $(HIST_OBJS) $(LIBS_PGPLOT)

mcmcphys: mcmcphys.f
	$(FC) $(FFLAGS) -o $(BIN)$@ mcmcphys.f

mcmc_corr: mcmc_corr.f
	$(FC) $(FFLAGS) -o $(BIN)$@ mcmc_corr.f $(LIBS_PGPLOT)

mcmchist5: mcmchist5.f $(HIST_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ mcmchist5.f $(HIST_OBJS) $(LIBS_PGPLOT)

HISTBF_OBJS = getrhostar.o histogram2.o avevar.o rqsort.o getfitpars.o ran2.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o
mcmchist5bf: mcmchist5bf.f $(HISTBF_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ mcmchist5bf.f $(HISTBF_OBJS) $(LIBS_PGPLOT)

# Object compilation rules
precision.o: $(UTILS)precision.f90
	$(FC) -c $(FFLAGS) -o $@ $<

fittermod.o: $(UTILS)fittermod.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

fittransitmodel3.o: $(UTILS)fittransitmodel3.f90 precision.o fittermod.o
	$(FC) -c $(FFLAGS) -o $@ $<

occultnlbm.o: $(UTILS)occultnlbm.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

readdata2.o: $(UTILS)readdata2.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

cutoutliers.o: $(UTILS)cutoutliers.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

meddiff.o: $(UTILS)meddiff.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

calcldprior.o: $(UTILS)calcldprior.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

%.o: $(UTILS)%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%.o: $(UTILS)%.f
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod *__genmod.f90 $(BIN)transitfit51 $(BIN)transitfit5 \
	      $(BIN)transitmcmc5 $(BIN)transitcut5 $(BIN)sigclip $(BIN)claretquad \
	      $(BIN)limbpriors $(BIN)dataclip $(BIN)transitsn5 $(BIN)transittiming5 \
	      $(BIN)prdump5 $(BIN)mcmcsetup $(BIN)transitdepth5 $(BIN)occultfind \
	      $(BIN)n1update $(BIN)datadump5 $(BIN)transitdur5 $(BIN)transitremove5 \
	      $(BIN)transitadd5 $(BIN)transitchisq5 $(BIN)transitsigclip5 \
	      $(BIN)transitscan5 $(BIN)rvplot $(BIN)plottimeseries $(BIN)starplot5 \
	      $(BIN)starplot5b $(BIN)transitplot5 $(BIN)updaten $(BIN)a2n \
	      $(BIN)mcmchistnew $(BIN)mcmcphys $(BIN)mcmc_corr $(BIN)mcmchist5 \
	      $(BIN)mcmchist5bf
