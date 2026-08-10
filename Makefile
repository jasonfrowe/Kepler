# ==============================================================================
# Top-level Makefile for Kepler Suite (Generated via ./configure or direct make)
# ==============================================================================

FC ?= gfortran
CC ?= gcc
OPT ?= -O3
OPENMP ?= -fopenmp
FFLAGS ?= $(OPT) $(OPENMP) -std=legacy -fallow-argument-mismatch
CFLAGS ?= $(OPT) -std=gnu89 -Wno-implicit-function-declaration

# Centralized bin directory
BIN ?= $(CURDIR)/bin/

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

# Subdirectories currently active
SUBDIRS = transitfit5 datatest transitfind

export FC CC OPT OPENMP FFLAGS CFLAGS BIN PGPLOT_DIR PGPLOT_LFLAGS X11_LFLAGS LIBS_PGPLOT LIBS_CFITSIO

.PHONY: all $(SUBDIRS) clean init_bin help

all: init_bin $(SUBDIRS)
	@echo "======================================================="
	@echo " All Kepler components successfully built into: $(BIN)"
	@echo "======================================================="

init_bin:
	@mkdir -p $(BIN)

transitfit5: init_bin
	$(MAKE) -C transitfit5 BIN=$(BIN)

datatest: init_bin
	$(MAKE) -C datatest BIN=$(BIN)

transitfind: init_bin
	$(MAKE) -C transitfind BIN=$(BIN)

clean:
	@for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean BIN=$(BIN); \
	done
	rm -rf $(BIN) *.o *.mod *__genmod.f90

help:
	@echo "Kepler Build System"
	@echo "Targets:"
	@echo "  all         - Build transitfit5, datatest, and transitfind (default)"
	@echo "  transitfit5 - Build only transitfit5 tools"
	@echo "  datatest    - Build only datatest tools"
	@echo "  transitfind - Build only transitfind tools"
	@echo "  clean       - Remove all compiled objects and binaries"
