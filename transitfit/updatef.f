      program updatef
C     Jason Rowe (2010) jasonfrowe@gmail.com
C     Used to update fitting parameters easily in batches
      implicit none
      integer nunit,nfit,iargc
      parameter(nfit=18)
      double precision sol(nfit),serr(nfit,2),Dpvary(nfit),err(nfit),
     .  doe,toff
      character*80 inputsol
      character*3 titles(nfit)

      if(iargc().lt.1) goto 901 !check number of command line arguements

      call getarg(1,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
      close(nunit) !release unit number as we are done with file
      
      serr(16,2)=-1.0d0 !turn on eclipse depth fitting
      
      write(0,*) "Exporting fit"     
      call exportfit(nfit,sol,serr,Dpvary,err,doe,toff,titles)
      
      
      goto 999
 901  write(0,*) "Usage: updatef f1.dat"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 999  end