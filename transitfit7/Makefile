#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F90 = ifort
#F90 = gfortran
#compiling object file flags
OPT1 = -O3 -fast -parallel -ipo
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
X11DIR = /usr/X11R6/lib
# libraries for linking PGPLOT
PLPLOTDIR = -I/usr/local/Cellar/plplot/5.9.11/lib/fortran/modules/plplot -I/usr/local/Cellar/plplot/5.9.11/include/plplot -L/usr/local/Cellar/plplot/5.9.11/lib 
LIBS = $(PLPLOTDIR) -lplplotf95d -lplplotf95cd
#LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng 
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/nbody/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: transitfit7 n2nb mcmcsetup7 transitmcmc7

MCMCINCL = precision.o readkeplc.o calcnbodies.o readinputsol.o ran2.o lcmodel.o opkdmain.o opkda1.o opkda2.o funcs.o transitmodel.o occultquad.o mandelagol.o dmcmc.o ovrwrt.o exportfit.o
transitmcmc7: transitmcmc7.f90 $(MCMCINCL)
	$(F90) $(LFLAGS) -o $(BIN)$@ transitmcmc7.f90 $(MCMCINCL)

MCMCSETUPINCL = precision.o fittingmod.o readkeplc.o calcnbodies.o readinputsol.o dfridr.o lcmodel.o funcs.o transitmodel.o occultquad.o mandelagol.o exportfit.o mercury.o
mcmcsetup7: mcmcsetup7.f90 $(MCMCSETUPINCL)
	$(F90) $(LFLAGS) -o $(BIN)$@ mcmcsetup7.f90 $(MCMCSETUPINCL)

TRANSINCL7 = precision.o fittingmod.o ocmod.o readkeplc.o calcnbodies.o readinputsol.o lcmodel.o aei2xyz.o rqsort.o funcs.o transitmodel.o occultquad.o mandelagol.o fittransitmodel.o minpack.o exportfit.o mercury.o octiming.o
transitfit7: transitfit7.f90 $(TRANSINCL7)
	$(F90) $(LFLAGS) -o $(BIN)$@ \
	transitfit7.f90 $(TRANSINCL7) 

N2NBINCL = precision.o readkeplc.o calcnbodies.o exportfit.o
n2nb: n2nb.f90 $(N2NBINCL)
	$(F90) $(LFLAGS) -o $(BIN)$@ n2nb.f90 $(N2NBINCL)

#building object libraries
ocmod.o: $(UTILS)ocmod.f90
	$(F90) $(FFLAGS) $(UTILS)ocmod.f90
octiming.o: $(UTILS)octiming.f90
	$(F90) $(FFLAGS) $(UTILS)octiming.f90
mercury.o: $(UTILS)mercury.f
	$(F90) $(FFLAGS) $(UTILS)mercury.f
rqsort.o: $(UTILS)rqsort.f
	$(F90) $(FFLAGS) $(UTILS)rqsort.f
aei2xyz.o: $(UTILS)aei2xyz.f90
	$(F90) $(FFLAGS) $(UTILS)aei2xyz.f90
ovrwrt.o: $(UTILS)ovrwrt.f
	$(F90) $(FFLAGS) $(UTILS)ovrwrt.f
dmcmc.o: $(UTILS)dmcmc.f90
	$(F90) $(FFLAGS) $(UTILS)dmcmc.f90
ran2.o: $(UTILS)ran2.f
	$(F90) $(FFLAGS) $(UTILS)ran2.f
dfridr.o: $(UTILS)dfridr.f
	$(F90) $(FFLAGS) $(UTILS)dfridr.f
fittingmod.o: $(UTILS)fittingmod.f90
	$(F90) $(FFLAGS) $(UTILS)fittingmod.f90
fittransitmodel.o: $(UTILS)fittransitmodel.f90
	$(F90) $(FFLAGS) $(UTILS)fittransitmodel.f90
exportfit.o: $(UTILS)exportfit.f90
	$(F90) $(FFLAGS) $(UTILS)exportfit.f90
mandelagol.o: $(UTILS)mandelagol.f
	$(F90) $(FFLAGS) $(UTILS)mandelagol.f
occultquad.o: $(UTILS)occultquad.f
	$(F90) $(FFLAGS) $(UTILS)occultquad.f
transitmodel.o: $(UTILS)transitmodel.f90
	$(F90) $(FFLAGS) $(UTILS)transitmodel.f90
lcmodel.o: $(UTILS)lcmodel.f90
	$(F90) $(FFLAGS) $(UTILS)lcmodel.f90
funcs.o: $(UTILS)funcs.f90
	$(F90) $(FFLAGS) $(UTILS)funcs.f90
readkeplc.o: $(UTILS)readkeplc.f90
	$(F90) $(FFLAGS) $(UTILS)readkeplc.f90
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
calcnbodies.o: $(UTILS)calcnbodies.f90
	$(F90) $(FFLAGS) $(UTILS)calcnbodies.f90
readinputsol.o: $(UTILS)readinputsol.f90
	$(F90) $(FFLAGS) $(UTILS)readinputsol.f90
opkdmain.o: $(UTILS)opkdmain.f
	$(F90) $(FFLAGS) $(UTILS)opkdmain.f
opkda1.o: $(UTILS)opkda1.f
	$(F90) $(FFLAGS) $(UTILS)opkda1.f
opkda2.o: $(UTILS)opkda2.f
	$(F90) $(FFLAGS) $(UTILS)opkda2.f
minpack.o: $(UTILS)minpack.f
	$(F90) $(FFLAGS) $(UTILS)minpack.f
	
# Removing object files
.PHONY : clean
clean :
	rm *.o
	rm *.mod
