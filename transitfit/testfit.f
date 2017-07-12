      program testfit
C     stupid program to debug that fucking retarded Mandel Agol program
      implicit none
      integer iargc,nunit,nmax,npt,nfit,i
      parameter(nmax=600000,nfit=16)
      double precision time(nmax),mag(nmax),merr(nmax),itime(nmax),
     .  Zerotime,sol(nfit),tmodel(nmax)
      character*80 obsfile
      data sol /1.0d0,1.0d0,1.0d0,1.0d0,3.52d0,88.0d0,0.0d0,0.0d0,
     .  0.0d0,0.4108d0,-0.1089d0,0.9040d0,-0.4374d0,0.0d0,0.0d0,0.0d0/
      

      if(iargc().lt.1) goto 901 !check number of command line arguements

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,mag,merr,itime,Zerotime)
      call readkeplc(nunit,nmax,npt,time,mag,merr,itime,Zerotime)
      close(nunit)!release unit number as we are done with the file.
      
      call transitmodel(npt,time,itime,tmodel,nfit,sol)
      write(0,*) (tmodel(i),i=1,npt)
      
      goto 999
 901  write(0,*) "Usage: transitfit datafile.dat"
      write(0,*) "  datafile.dat : contains MOST photometry"
      goto 999
 903  write(0,*) "Cannot open ",obsfile
      goto 999
 999  end