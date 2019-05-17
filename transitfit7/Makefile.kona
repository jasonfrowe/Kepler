#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F90 = ifort
F77 = ifort
#F90 = gfortran
#compiling object file flags
OPT1 = -O3
#OPT1 = -O0 -g -CB
#OPT1 = -O3 
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
X11DIR = /usr/X11/lib
# libraries for linking PGPLOT
PLPLOTDIR = -I/usr/local/Cellar/plplot/5.9.11/lib/fortran/modules/plplot -I/usr/local/Cellar/plplot/5.9.11/include/plplot -L/usr/local/Cellar/plplot/5.9.11/lib 
LIBS = $(PLPLOTDIR) -lplplotf95d -lplplotf95cd
#LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng 
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/nbody/bin/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: transitfit7 n2nb mcmcsetup7 transitmcmc7

MCMCINCL = precision.o readkeplc.o calcnbodies.o readinputsol.o ran2.o lcmodel.o opkdmain.o \
opkda1.o opkda2.o funcs.o transitmodel.o occultquad.o mandelagol.o dmcmc.o ovrwrt.o exportfit.o \
aei2xyz.o rqsort.o mercury.o
transitmcmc7: transitmcmc7.f90 $(MCMCINCL)
	$(F90) $(LFLAGS) -o $(BIN)$@ transitmcmc7.f90 $(MCMCINCL)

MCMCSETUPINCL = precision.o fittingmod.o readkeplc.o calcnbodies.o readinputsol.o dfridr.o lcmodel.o funcs.o transitmodel.o \
occultquad.o mandelagol.o exportfit.o mercury.o aei2xyz.o rqsort.o ran2.o
mcmcsetup7: mcmcsetup7.f90 $(MCMCSETUPINCL)
	$(F90) $(LFLAGS) -o $(BIN)$@ mcmcsetup7.f90 $(MCMCSETUPINCL)

TRANSINCL7 = precision.o fittingmod.o ocmod.o readkeplc.o calcnbodies.o readinputsol.o lcmodel.o \
aei2xyz.o rqsort.o transitmodel.o occultquad.o mandelagol.o minpack.o exportfit.o keplerian.o \
octiming.o ocmod.o percorcalc.o mercury.o median.o fittransitmodel.o
#TRANSINCL7 = precision.o fittingmod.o ocmod.o readkeplc.o calcnbodies.o readinputsol.o lcmodel.bs_v1.o \
#aei2xyz2.o rqsort.o transitmodel.o occultquad.o mandelagol.o minpack.o exportfit.o keplerian.o \
#octiming.o ocmod.o percorcalc.o median.o bsstep.o bsstep_f.o
transitfit7: transitfit7.f90 $(TRANSINCL7)
	$(F90) $(LFLAGS) -o $(BIN)$@ \
	transitfit7.f90 $(TRANSINCL7) $(LIBS2)

N2NBINCL = precision.o readkeplc.o calcnbodies.o exportfit.o
n2nb: n2nb.f90 $(N2NBINCL)
	$(F90) $(LFLAGS) -o $(BIN)$@ n2nb.f90 $(N2NBINCL)

#building object libraries
#building object libraries
%_mp.o : $(UTILS)%.f90
	$(F90) $(FFLAGS) $(LMP) -o $@ $(UTILS)$*.f90

%_mp.o : $(UTILS)%.f
	$(F77) $(FFLAGS) $(LMP) -o $@ $(UTILS)$*.f

%.o : $(UTILS)%.f90
	$(F90) $(FFLAGS) -o $@ $<

%.o : $(UTILS)%.f
	$(F77) $(FFLAGS) -o $@ $<

# Removing object files
.PHONY : clean
clean :
	rm *.o
	rm *.mod
