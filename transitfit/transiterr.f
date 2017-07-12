      program transiterr
C     determines uncertainities from transitfits
      integer iargc,nunit,nmax,npt,nfit,npar
      parameter(nmax=650000,nfit=16)
      double precision time(nmax),mag(nmax),merr(nmax),itime(nmax),
     .  Keplertime,sol(nfit),serr(nfit,2),Dpvary(nfit),toff,sol2,
     .  tmodel(nmax)
      character*80 obsfile,inputsol



      if(iargc().lt.2) goto 901 !check number of command line arguements

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,mag,merr,itime,Keplertime)
      call readkeplc(nunit,nmax,npt,time,mag,merr,itime,Keplertime)
      close(nunit)!release unit number as we are done with the file.

      call getarg(2,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,sol,serr,Dpvary,toff)
      close(nunit) !release unit number as we are done with file
      
      npar=1
      call marginalize(npt,time,mag,merr,itime,tmodel,nfit,sol,
     .  serr,Dpvary,npar,sol2)
      
      goto 999
 901  write(0,*) "Usage: transiterr datafile.dat fitpars.dat"
      write(0,*) "  datafile.dat : contains MOST photometry"
      write(0,*) "  fitspars.dat : contains fitting information"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 903  write(0,*) "Cannot open ",obsfile
 999  end