      program transitremove
      implicit none
      integer iargc,nunit,nmax,npt,nfit,i
      parameter(nmax=600000,nfit=18)
      integer dtype(nmax)
      double precision time(nmax),mag(nmax),merr(nmax),itime(nmax),
     .  MOSTtime,sol(nfit),serr(nfit,2),err(nfit),doe,toff,Dpvary(nmax),
     .  tmodel(nmax)
      character*80 obsfile,inputsol


      if(iargc().lt.2) goto 901 !check number of command line arguements

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      call readkeplc(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      close(nunit)!release unit number as we are done with the file.

      call getarg(2,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
      close(nunit) !release unit number as we are done with file
      
      do 17 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
 17   continue     
      
      call transitmodel(npt,time,itime,dtype,tmodel,nfit,sol)
      
      do 20 i=1,npt
        write(6,*) time(i)+54900.0d0-0.5d0,
     .          10.0**((mag(i)-tmodel(i))/-2.5d0)-1.0d0,merr(i)
c        write(6,*) time(i)+54900.0d0-0.5d0,
c     .          10.0**((mag(i)+tmodel(i))/-2.5d0)-1.0d0,merr(i)
 20   continue
      
      goto 999
 901  write(0,*) "Usage: transitremove photfile f1.dat"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 903  write(0,*) "Cannot open ",obsfile
      goto 999
 999  end