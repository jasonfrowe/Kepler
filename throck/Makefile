#Name of C compiler
CCMPL = icc
#Name of Fortran compiler
F77 = ifort
#compiling object file flags
OPT1 = -O3 
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
X11DIR = /usr/X11/lib
# libraries for linking PGPLOT
LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/bin/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: onedrock twodrock threedrock

THROCKDIR = $(BIN)
onedrock: onedrock.f 
	$(F77) $(LFLAGS) -o $(THROCKDIR)$@ onedrock.f $(LIBS)

TWODINCL = surfacearea.o keplerian.o stellarrad.o thplot.o radiate.o conduction.o heatcapacity.o
twodrock: twodrock.f $(TWODINCL)
	$(F77) $(LFLAGS) -o $(THROCKDIR)$@ twodrock.f \
	$(TWODINCL) $(LIBS)

THREEDINCL = surfacearea.o thplot.o keplerian.o stellarrad3D.o heatcapacity3D.o radiate3D.o conduction3D.o kcp3D.o keplerbandpass.o planck.o genfluxes.o spline.o obsflux.o getrockpars.o
threedrock: threedrock.f $(THREEDINCL)
	$(F77) $(LFLAGS) -o $(THROCKDIR)$@ threedrock.f \
	$(THREEDINCL) $(LIBS)

#building object libraries

keplerian.o: common/keplerian.f
	$(F77) $(FFLAGS) common/keplerian.f
surfacearea.o: common/surfacearea.f
	$(F77) $(FFLAGS) common/surfacearea.f
stellarrad.o: 2D/stellarrad.f
	$(F77) $(FFLAGS) 2D/stellarrad.f
thplot.o: common/thplot.f
	$(F77) $(FFLAGS) common/thplot.f
radiate.o: 2D/radiate.f
	$(F77) $(FFLAGS) 2D/radiate.f
conduction.o: 2D/conduction.f
	$(F77) $(FFLAGS) 2D/conduction.f
heatcapacity.o: 2D/heatcapacity.f
	$(F77) $(FFLAGS) 2D/heatcapacity.f
stellarrad3D.o: 3D/stellarrad3D.f
	$(F77) $(FFLAGS) 3D/stellarrad3D.f
heatcapacity3D.o: 3D/heatcapacity3D.f
	$(F77) $(FFLAGS) 3D/heatcapacity3D.f
radiate3D.o: 3D/radiate3D.f
	$(F77) $(FFLAGS) 3D/radiate3D.f
conduction3D.o: 3D/conduction3D.f
	$(F77) $(FFLAGS) 3D/conduction3D.f
kcp3D.o: 3D/kcp3D.f
	$(F77) $(FFLAGS) 3D/kcp3D.f
keplerbandpass.o: flux/keplerbandpass.f
	$(F77) $(FFLAGS) flux/keplerbandpass.f
planck.o: flux/planck.f
	$(F77) $(FFLAGS) flux/planck.f
genfluxes.o: flux/genfluxes.f
	$(F77) $(FFLAGS) flux/genfluxes.f
spline.o: common/spline.f
	$(F77) $(FFLAGS) common/spline.f
obsflux.o: 3D/obsflux.f
	$(F77) $(FFLAGS) 3D/obsflux.f
getrockpars.o: 3D/getrockpars.f
	$(F77) $(FFLAGS) 3D/getrockpars.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
