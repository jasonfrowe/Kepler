CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getpars(nunit,parfile,dtime,tmax,Re,PorbD,PerD,asemiAU,
     .  eccn,obliquity,alphad,Lstar,Rstar,Ab,rdepth,rho,tdr,Tstar,Eint)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit
      double precision dtime,tmax,Re,PorbD,PerD,asemiAU,eccn,obliquity,
     .  alphad,Lstar,Rstar,Ab,rdepth,rho,tdr,Tstar,Eint
      character*80 parfile
      
      open(unit=nunit,file=parfile,status='old',err=901)
      read(nunit,*) dtime
      read(nunit,*) tmax
      read(nunit,*) Re
      read(nunit,*) PorbD
      read(nunit,*) PerD
      read(nunit,*) asemiAU
      read(nunit,*) eccn
      read(nunit,*) obliquity
      read(nunit,*) alphad
      read(nunit,*) Lstar
      read(nunit,*) Rstar
      read(nunit,*) Ab
      read(nunit,*) rdepth
      read(nunit,*) rho
      read(nunit,*) tdr
      read(nunit,*) Tstar
      read(nunit,*) Eint
      
      close(nunit)
      goto 999

C     Default Model parameters 
c      dtime=30.0d0 !set time steps
c      tmax=100.0*365.25*24.0*60.0*60.0 !stop model at this time
c      Re=1.0d0 !radius of rock in Earth Radii
c      PorbD=365.25d0 !orbital period of planet in days
c      PerD=1.0d0 !rotational period of planet in days
c      asemiAU=1.0d0 !semi-major axis [AU]
c      eccn=0.016710219d0 !eccentricity 
c      obliquity=23.439281d0 !obliquity of planet rotation
c      alphad=90.0d0 !angle of rotated grid [degrees] - leave alone
c      Lstar=1.0d0 !Luminosity of star [Lsun]
c      Rstar=1.0d0 !Radius of star [Rsun]
c      Ab=0.367 !Bond Albedo of planet surface
c      rdepth=10.0d0 !depth of layer [m]
c      rho=2700.0d0 !density [kg m-3]
c      tdr=20.0d0 !Temperature change / unit depth [K km-1]
c      Tstar=5780.0d0 !star temperature [K]
c      Eint=0.078d0 !Internal Energy output of planet [W m-2]
      
 901  write(0,*) "Cannot open ",parfile
      pause     
      goto 999
 999  return
      end