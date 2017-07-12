      program transitduration
      implicit none
C     Calculate transit duration from Transit model
      integer nunit,nfit
      parameter(nfit=18)
      double precision sol(nfit),serr(nfit,2),Dpvary(nfit),err(nfit),
     .  doe,toff,Pi,tPi,transitdur,tdur
      character*80 inputsol

      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi   
      
      if(iargc().lt.1) goto 901 !check number of command line arguements

      call getarg(1,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
      close(nunit) !release unit number as we are done with file
      
      tdur=transitdur(nfit,sol)
      write(6,500) tdur/3600.0d0
 500  format(F7.4)     
      
      goto 999
 901  write(0,*) "Usage: transitdur f1.dat"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdur(nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This is for circular orbits only.
      implicit none
      integer nfit
      double precision sol(nfit),Psec,M1,M2,R1,R2,asemi,temp(4),incl
C     Read in physical constants (Pi,G,Msun,..etc.)      
      include "utils/physcons.f"
      
      M1=sol(1)*Msun !kg ; mass of star
      M2=sol(2)*Mjup !kg ; mass of planet
      Psec=sol(5)*24.0*60.0*60.0 !sec ; period of planet
      asemi=(Psec*Psec*G*(M1+M2)/(4.0*Pi*Pi))**(1.0/3.0) !m
      R1=sol(3)*Rsun  !radius of star
      R2=sol(4)*Rjup  !radius of planet
      incl=Pi*sol(6)/180.0d0
      
      temp(1)=Psec/Pi
      temp(2)=R1/asemi
      temp(3)=(1+(R2/R1))**2.0-((asemi/R1)*cos(incl))**2.0
      temp(4)=1-cos(incl)*cos(incl)     
      
      transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))
      
      return
      end