      subroutine getinputpars(nunit,nfrho,adrsi,adrsierr,massi,radi,
     .  agei,Teff,Tefferr,Z,Zerr,dZ,L,Lerr,Psec,dmass,dage,workdir,rhoc,
     .  rhocerr,logg,loggerr,aotu)
      implicit none
      integer nunit, ndum,nfrho
      double precision adrsi,adrsierr(6),massi,radi,Teff,Tefferr,rdum,
     .  Z,Zerr(2),dmass,dage,Per,Psec,rhoc,rhocerr,dZ,adrssdev,L,Lerr,
     .  logg,loggerr,G,Msun,Rsun,rhoi,Pi,FeH,FeHerr,Zsun,dum,aotu,agei
      character*80 workdir
      
      Zsun=0.01812
      
      nfrho=0 !use adrs to fit rho_*
      read(nunit,*) ndum,rdum,adrsi,adrssdev,adrsierr(4),adrsierr(3),
     .  adrsierr(5),adrsierr(2),adrsierr(6),adrsierr(1)
      if(adrssdev.le.0.0d0) nfrho=1 !set flag to not fit rho_*
      read(nunit,*) Teff,Tefferr
c      read(nunit,*) Z,Zerr,dZ !metals, error in measurement, mcmc tuning
      read(nunit,*) FeH,FeHerr,dZ
      read(nunit,*) L,Lerr !stellar luminosity, error
      read(nunit,*) logg,loggerr !constrains on log(g)
      read(nunit,*) Per
      Psec=per*24.0d0*60.0d0*60.0d0 !period (sec) 
c      read(nunit,*) massi !initial mass
c      read(nunit,*) radi !initial radius
      read(nunit,*) dmass !mcmc tuning for mass
      read(nunit,*) dage  !mcmc tuning for age
      read(nunit,*) rhoc,rhocerr !companion density
      read(nunit,*) workdir  !completely useless.
      
c      Z=Zsun*10**(FeH)
c      Zerr(1)=Zsun*10**(FeH+FeHerr)-Zsun*10**(FeH)
c      Zerr(2)=Zsun*10**(FeH)-Zsun*10**(FeH-FeHerr)

      call feh2z(FeH,Z)
      dum=FeH+FeHerr
      call feh2z(dum,Zerr(1))
      Zerr(1)=Zerr(1)-Z
      dum=FeH-FeHerr
      call feh2z(dum,Zerr(2))
      Zerr(2)=Z-Zerr(2)
      
      if(FeHerr.eq.0.0d0) then
        Zerr(1)=0.0d0
        dZ=0.0
      endif
      write(0,*) "Z:   ",Z
      write(0,*) "Zerr ",Zerr(1),Zerr(2)
            
      G=6.67259E-11
      Rsun=6.9599E8
      Msun=1.989E30
      Pi=acos(-1.d0)!define Pi and 2*Pi
      write(0,*) "Start mass scan"
      call massscan(massi,radi,agei,Teff,Tefferr,logg,loggerr,Psec,
     .  L,Lerr,adrsi,adrsierr,rhoc,rhocerr,Z,nfrho,aotu)
      write(0,*) "Stop mass scan"

c      call msmasstemp(Teff,massi)
c      write(0,*) "massi",massi
c      if(nfrho.eq.0)then
c        rhoi=adrsi**3.0*Pi*3.0d0/(Psec*Psec*G)-rhoc
c        radi= (massi*Msun/(4.0/3.0*pi*rhoi))**(1.0/3.0)/Rsun
c      elseif(loggerr.gt.0.0d0)then
cc        write(0,*) massi,sqrt(G*massi*Msun/10**(logg-2.0d0)),Rsun
c        radi= sqrt(G*massi*Msun/10**(logg-2.0d0))/Rsun
c      else
c        radi=massi**0.8 !M-S Mass-radius relation
c      endif
      write(0,*) "mass,radi:",massi,radi
c      read(5,*)
      
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine feh2z(FeH,Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (dydz=2.0E0)
      parameter (yp=0.23E0)
      parameter (zp=0.0E0)
      parameter (Zsun=0.0181D0)
      parameter (Xsun=0.7148705D0)
      parameter (FeHa2=-0.217D0)
      parameter (FeHa4=-0.470D0)
      common /head/YYalpfe

      quad(x1,y1,x2,y2,x3,y3,x)=y1*(x2-x)*(x3-x)/((x2-x1)*(x3-x1))
     +                         +y2*(x1-x)*(x3-x)/((x1-x2)*(x3-x2))
     +                         +y3*(x1-x)*(x2-x)/((x1-x3)*(x2-x3))

       FeH0=FeH
     +      -quad(0.0d0,0.0d0,0.3d0,FeHa2,0.6d0,FeHa4,YYalpfe)
C       if(abs(YYalpfe-0.3d0).le.0.1d0)FeH0=FeH-FeHa2
C       if(abs(YYalpfe-0.6d0).le.0.1d0)FeH0=FeH-FeHa4
      ZovX=(10.0d0**FeH0)*Zsun/Xsun
      Z=ZovX*(1.0d0+dydz*zp-yp)/(1.0d0+ZovX*(1.0d0+dydz))
      return
      end