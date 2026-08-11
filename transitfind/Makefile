# ==============================================================================
# Makefile for transitfind
# ==============================================================================

FC ?= gfortran
CC ?= gcc
OPT ?= -O3
OPENMP ?= -fopenmp

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
  LDFLAGS ?= -Wl,-stack_size,0x20000000
  X11_LFLAGS ?= -L/opt/homebrew/lib -L/opt/X11/lib -lX11 -lpng -lz
else
  X11_LFLAGS ?= -L/usr/X11R6/lib -L/usr/lib64 -L/usr/local/lib -lX11 -lpng -lz
endif

FFLAGS ?= $(OPT) $(OPENMP) $(LDFLAGS) -std=legacy -fallow-argument-mismatch
CFLAGS ?= $(OPT) -std=gnu89 -Wno-implicit-function-declaration

LIBS_PGPLOT = $(PGPLOT_LFLAGS) $(X11_LFLAGS)

# Targets to build
TARGETS = transitfind5 transitfind3 transitfind2 transitfind pdmsearch \
          pickcandidates tpsclean tpsdawg tpsstats effplot notesmerge mescomp

.PHONY: all clean init_bin

all: init_bin $(TARGETS)

init_bin:
	@mkdir -p $(BIN)

# Common object dependencies
TF_OBJS = precision.o rqsort.o stdev.o graphics.o

transitfind5: transitfind5.f90 $(TF_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ transitfind5.f90 $(TF_OBJS) $(LIBS_PGPLOT)

transitfind3: transitfind3.f90 $(TF_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ transitfind3.f90 $(TF_OBJS) $(LIBS_PGPLOT)

transitfind2: transitfind2.f90 $(TF_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ transitfind2.f90 $(TF_OBJS) $(LIBS_PGPLOT)

TF1_OBJS = dataio.o bls.o nyquest.o rqsort.o
transitfind: transitfind1.f $(TF1_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ transitfind1.f $(TF1_OBJS) $(LIBS_PGPLOT)

PDM_OBJS = precision.o readdata2.o nyquest2.o rqsort.o calcnsteps.o qfunc.o \
           pdm.o ovrwrt.o makeplot.o graphics.o
pdmsearch: pdmsearch.f90 $(PDM_OBJS)
	$(FC) $(FFLAGS) -o $(BIN)$@ pdmsearch.f90 $(PDM_OBJS) $(LIBS_PGPLOT)

pickcandidates: pickcandidates.f90 precision.o
	$(FC) $(FFLAGS) -o $(BIN)$@ pickcandidates.f90 precision.o

tpsclean: tpsclean.f
	$(FC) $(FFLAGS) -o $(BIN)$@ tpsclean.f

tpsdawg: tpsdawg.f
	$(FC) $(FFLAGS) -o $(BIN)$@ tpsdawg.f $(LIBS_PGPLOT)

tpsstats: tpsstats.f
	$(FC) $(FFLAGS) -o $(BIN)$@ tpsstats.f $(LIBS_PGPLOT)

effplot: effplot.f
	$(FC) $(FFLAGS) -o $(BIN)$@ effplot.f $(LIBS_PGPLOT)

notesmerge: notesmerge.f
	$(FC) $(FFLAGS) -o $(BIN)$@ notesmerge.f

mescomp: mescomp.f
	$(FC) $(FFLAGS) -o $(BIN)$@ mescomp.f $(LIBS_PGPLOT)

# Object compilation rules
precision.o: $(UTILS)precision.f90
	$(FC) -c $(FFLAGS) -o $@ $<

readdata2.o: $(UTILS)readdata2.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

nyquest2.o: $(UTILS)nyquest2.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

calcnsteps.o: $(UTILS)calcnsteps.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

qfunc.o: $(UTILS)qfunc.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

pdm.o: $(UTILS)pdm.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

makeplot.o: $(UTILS)makeplot.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

makeplot2.o: $(UTILS)makeplot2.f90 precision.o
	$(FC) -c $(FFLAGS) -o $@ $<

%.o: $(UTILS)%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%.o: $(UTILS)%.f
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod *__genmod.f90 $(BIN)transitfind5 $(BIN)transitfind3 \
	      $(BIN)transitfind2 $(BIN)transitfind $(BIN)pdmsearch \
	      $(BIN)pickcandidates $(BIN)tpsclean $(BIN)tpsdawg $(BIN)tpsstats \
	      $(BIN)effplot $(BIN)notesmerge $(BIN)mescomp
