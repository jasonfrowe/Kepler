#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
FPATH = /opt/intel/compilers_and_libraries_2017.2.163/mac/bin/intel64/
F77 = $(FPATH)ifort
F90 = $(FPATH)ifort
#F77 = gfortran
#compiling object file flags
OPT1 = -O3 -parallel -ipo 
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
X11DIR = /usr/X11/lib
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

all: tempplot teqplot plot975 teffrhostar

hzfigure1incl = precision.o getfitpars.o transitmodel.o mandelagol.o keplerian.o occultquad.o ttcor.o phasept.o binp.o sort2.o 
hzfigure1inclbin = /Users/rowe/Documents/Kepler/hzplanets/
hzfigure1: hzfigure1.f90 $(hzfigure1incl)
	$(F90) $(LFLAGS) -o $(hzfigure1inclbin)$@ hzfigure1.f90 $(hzfigure1incl) $(LIBS)

k10phaseplotincl = precision.o rqsort.o
k10phaseplotbin = /Users/rowe/Documents/transitfit/Kepler_n/koi72.n/
k10phaseplot: k10phaseplot.f90 $(k10phaseplotincl)
	$(F90) $(LFLAGS) -o $(k10phaseplotbin)$@ k10phaseplot.f90 $(k10phaseplotincl) $(LIBS)

teffrhostarincl = 
teffrhostarbin = /Users/rowe/Documents/Kepler/fpwg/
teffrhostar: teffrhostar.f90 $(teffrhostarincl)
	$(F90) $(LFLAGS) -o $(teffrhostarbin)$@ teffrhostar.f90 $(teffrhostarincl) $(LIBS)

koi3138incl = precision.o ran2.o
koi3138bin = /Users/rowe/Documents/transitfit/Kepler_n/koi3138.n/
koi3138srad: koi3138srad.f90 $(koi3138incl)
	$(F90) $(LFLAGS) -o $(koi3138bin)$@ koi3138srad.f90 $(koi3138incl) $(LIBS)

koi3138spec: koi3138spec.f90 $(koi3138incl)
	$(F90) $(LFLAGS) -o $(koi3138bin)$@ koi3138spec.f90 $(koi3138incl) $(LIBS)


KOI314INCL = precision.o
KOI314bin = /Users/rowe/Documents/transitfit/Kepler_n/koi314.n/
koi314rho: koi314rho.f90 $(KOI314INCL)
	$(F90) $(LFLAGS) -o $(KOI314bin)$@  koi314rho.f90 $(KOI314INCL) $(LIBS)

koi81_track: koi81_track.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/Kepler_n.old/koi81.n/tracks/$@ koi81_track.f $(LIBS)
koi701_track: koi701_track.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/Kepler_n/koi701.n/tracks/$@ koi701_track.f $(LIBS)
rockyvsgasfig: rockyvsgasfig
	$(F77) $(LFLAGS) -o /home/rowe/Documents/tex2/rocksvsgas/$@ rockyvsgasfig.f $(LIBS)
tempalbedo: tempalbedo.f
	$(F77) $(LFLAGS) -o $(BIN)$@ tempalbedo.f $(LIBS)
teqplot2: teqplot2.f
	$(F77) $(LFLAGS) -o /home/rowe/Kepler/Teq/$@ teqplot2.f $(LIBS)
kicupdateplot: kicupdateplot.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/smodels/YY/$@ kicupdateplot.f $(LIBS)
koi70bcdplot: koi70bcdplot.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/Kepler_new/koi70.n/$@ koi70bcdplot.f $(LIBS)
koi70p4plot: koi70p4plot.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/Kepler_new/koi70.n/$@ koi70p4plot.f $(LIBS)
nfitplots: nfitplots.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/$@ nfitplots.f $(LIBS)
koi707_trackv2: koi707_trackv2.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/Kepler_new/koi707.n/tracks/$@ koi707_trackv2.f $(LIBS)
koi707_track: koi707_track.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/Kepler_new/koi707.n/tracks/$@ koi707_track.f $(LIBS)
koi70rvplot: koi70rvplot.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/Kepler_new/koi70.n/$@ koi70rvplot.f $(LIBS)
koi87_track: koi87_track.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/Kepler/koi87/tracks/$@ koi87_track.f $(LIBS)

koi1728_track: koi1728_track.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/Kepler/koi1728/tracks/$@ koi1728_track.f $(LIBS)

koi20_track: koi20_track.f
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/Kepler/koi20/tracks/$@ koi20_track.f $(LIBS)

PLOT975INCL = gettrack.o yytrack.o
plot975: plot975.f $(PLOT975INCL)
	$(F77) $(LFLAGS) -o /home/rowe/p555/transitfit/Kepler/koi975/plot/$@ plot975.f $(PLOT975INCL) $(LIBS)
tempplot: tempplot.f
	$(F77) $(LFLAGS) -o $(BIN)$@ tempplot.f $(LIBS)
teqplot: teqplot.f
	$(F77) $(LFLAGS) -o /home/rowe/Kepler/Teq/$@ teqplot.f $(LIBS)

#building object libraries

sort2.o: $(UTILS)sort2.f
	$(F77) $(FFLAGS) $(UTILS)sort2.f
binp.o: $(UTILS)binp.f
	$(F77) $(FFLAGS) $(UTILS)binp.f
phasept.o: $(UTILS)phasept.f
	$(F77) $(FFLAGS) $(UTILS)phasept.f
transitmodel.o: $(UTILS)transitmodel.f
	$(F77) $(FFLAGS) $(UTILS)transitmodel.f
mandelagol.o: $(UTILS)mandelagol.f
	$(F77) $(FFLAGS) $(UTILS)mandelagol.f
occultquad.o: $(UTILS)occultquad.f
	$(F77) $(FFLAGS) $(UTILS)occultquad.f
keplerian.o: $(UTILS)keplerian.f
	$(F77) $(FFLAGS) $(UTILS)keplerian.f
ttcor.o: $(UTILS)ttcor.f
	$(F77) $(FFLAGS) $(UTILS)ttcor.f
getfitpars.o: $(UTILS)getfitpars.f
	$(F77) $(FFLAGS) $(UTILS)getfitpars.f
rqsort.o: $(UTILS)rqsort.f
	$(F77) $(FFLAGS) $(UTILS)rqsort.f
ran2.o: $(UTILS)ran2.f
	$(F77) $(FFLAGS) $(UTILS)ran2.f
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
freqout.o: $(UTILS)freqout.f
	$(F77) $(FFLAGS) $(UTILS)freqout.f
gettrack.o: $(UTILS)gettrack.f
	$(F77) $(FFLAGS) $(UTILS)gettrack.f
yytrack.o: $(UTILS)yytrack.f
	$(F77) $(FFLAGS) $(UTILS)yytrack.f

# Removing object files
.PHONY : clean
clean :
	rm *.o
