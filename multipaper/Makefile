#Name of C compiler
CCMPL = icc
#Name of Fortran compiler
F77 = ifort
#F77 = gfortran
#compiling object file flags
OPT1 = -O3 -fast -parallel -ipo -no-prec-div
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
X11DIR = /usr/X11R6/lib
# libraries for linking PGPLOT
LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng 
#LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng 
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/Kepler/multi/
#utils source directory
UTILS = utils/

#Listing of programs to create.

all: targetlist mescomp multistats cumldisp multistatsV2 singles quarterobs tdurcomp radectable tdurplots eccndisp popplot histanim

histanim: histanim.f
	$(F77) $(LFLAGS) -o $(BIN)media/$@ histanim.f $(LIBS)
popplot: popplot.f
	$(F77) $(LFLAGS) -o $(BIN)$@ popplot.f $(LIBS)
eccndisp: eccndisp.f
	$(F77) $(LFLAGS) -o $(BIN)$@ eccndisp.f $(LIBS)
tdurplots: tdurplots.f
	$(F77) $(LFLAGS) -o $(BIN)$@ tdurplots.f $(LIBS)
radectable: radectable.f
	$(F77) $(LFLAGS) -o $(BIN)$@ radectable.f
tdurcomp: tdurcomp.f
	$(F77) $(LFLAGS) -o $(BIN)$@ tdurcomp.f $(LIBS)
quarterobs: quarterobs.f
	$(F77) $(LFLAGS) -o $(BIN)$@ quarterobs.f
singles: singles.f
	$(F77) $(LFLAGS) -o $(BIN)$@ singles.f $(LIBS)
multistatsV2: multistatsV2.f rqsort.o
	$(F77) $(LFLAGS) -o $(BIN)$@ multistatsV2.f $(LIBS) rqsort.o
cumldisp: cumldisp.f
	$(F77) $(LFLAGS) -o $(BIN)$@ cumldisp.f $(LIBS)

multistats: multistats.f
	$(F77) $(LFLAGS) -o $(BIN)$@ multistats.f $(LIBS)

mescomp: mescomp.f
	$(F77) $(LFLAGS) -o $(BIN)$@ mescomp.f $(LIBS)

targetlistincl = 
targetlist: targetlist.f $(targetlistincl)
	$(F77) $(LFLAGS) -o $(BIN)$@ targetlist.f $(targetlistincl)

#building object libraries

starplot.o: $(UTILS)starplot.f
	$(F77) $(FFLAGS) $(UTILS)starplot.f
	
# Removing object files
.PHONY : clean
clean :
	rm *.o
rqsort.o: $(UTILS)rqsort.f
	$(F77) $(FFLAGS) $(UTILS)rqsort.f
