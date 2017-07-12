      program stellardensity
      implicit none
      integer nunit,nfit,iargc
      parameter(nfit=18)
      double precision sol(nfit),serr(nfit,2),err(nfit),adrs,per,tpi,
     .  eccn,Pi,pid2,rhos,G,Psec
      character*80 inputsol
      
      if(iargc().lt.1) goto 901

      call getarg(1,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,sol,serr,err)
      close(nunit) !release unit number as we are done with file   

      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi 
      pid2=Pi/2.0d0
      G=6.674d-11 !N m^2 kg^-2  Gravitation constamt
      
      per=sol(2)     !Period (days)   
      Psec=per*24.0d0*60.0d0*60.0d0 !period (sec)   
      eccn=sqrt(sol(7)*sol(7)+sol(8)*sol(8)) !eccentricity
      if(eccn.ge.1.0) eccn=0.99 
C     a/R*
c      adrs=sol(5)*per/tpi*sqrt(1-sol(3))*(1+sol(8))/sqrt(1-eccn*eccn)
      adrs=sol(5)*per/tpi*sqrt((1.0d0+sol(4))**2.0d0-sol(3))*
     .  (1.0+sol(8))/sqrt(1.0-eccn*eccn)
      rhos= adrs**3.0*Pi*3.0d0/(Psec*Psec*G)
      write(6,*) "Stellar density: ",rhos/1000.0d0
      
      goto 999
 901  write(0,*) "Usage: stellardensity a1.dat"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 999  end