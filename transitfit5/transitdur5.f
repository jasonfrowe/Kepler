CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program transitdur5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Fits for mean stellar density (p_star) for multi-planet candidates
      implicit none
      integer iargc,nfit,i,nunit,nplanet,np
      parameter(nfit=108)
      double precision sol(nfit),serr(nfit,2),err(nfit),tdur,transitdur,
     .  transitdurflat,tdurf
      character*80 inputsol,cline
     
      if(iargc().lt.2) goto 901

      call getarg(1,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
c      write(0,*) "reading in input solution"
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
c      write(0,*) "done reading input solution"
      close(nunit) !release unit number as we are done with file
          
      call getarg(2,cline)
      read(cline,*) np
      if((np.lt.1).or.(np.gt.nplanet))goto 901
      
      tdur=transitdur(nfit,sol,np) !returns transit duration in hours
c      tdurf=transitdurflat(nfit,sol,np) !'flat' part of transits
      write(6,500) tdur
! 500  format(F7.4)
 500  format(F13.9)
       
      goto 999
 901  write(0,*) "Usage: transitdur5 <n1.dat> <np>"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdur(nfit,sol,np)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,col,np
      double precision sol(nfit),Pi,Msun,Rsun,G,aConst,b,Psec,adrs,
     .  cincl,temp(4),tpi,rdr,bb

      Pi=acos(-1.d0)   !Pi
      tpi=2.0d0*pi
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

      col=8+(np-1)*10
      b=sol(col+3)
      bb=b*b
      Psec=sol(col+2)*8.64d4 !sec ; period of planet
      adrs=1000.0*sol(1)*G*(Psec)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0) !a/R*
      cincl=b/adrs !cos(i)
      rdr=sol(col+4)
        
      temp(1)=Psec/Pi
      temp(2)=1.0d0/adrs
      temp(3)=(1+rdr)**2.0-bb
      temp(4)=1-cincl*cincl
C     Transit duration in hours
      transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))/3600.0
      
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdurflat(nfit,sol,np)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,col,np
      double precision sol(nfit),Pi,Msun,Rsun,G,aConst,b,Psec,adrs,
     .  cincl,temp(4),tpi,rdr,bb

      Pi=acos(-1.d0)   !Pi
      tpi=2.0d0*pi
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

      col=8+(np-1)*10
      b=sol(col+3)
      bb=b*b
      Psec=sol(col+2)*8.64d4 !sec ; period of planet
      adrs=1000.0*sol(1)*G*(Psec)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0) !a/R*
      cincl=b/adrs !cos(i)
      rdr=sol(col+4)
        
      temp(1)=Psec/Pi
      temp(2)=1.0d0/adrs
      temp(3)=(1-rdr)**2.0-bb
      temp(4)=1-cincl*cincl
C     Transit duration in hours
      transitdurflat=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))/3600.0
      
      
      return
      end
      
