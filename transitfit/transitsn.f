      program transitSN
      implicit none
C     Convert a Physical space model to a geometric model
      integer nfitg,nfit,nunit,iargc,i,npt,nmax,npt1,j
      parameter(nfitg=18,nfit=18,nmax=600000)
      double precision solg(nfitg),serrg(nfitg,2),Dpvaryg(nfitg),
     .  errg(nfitg),doe,toff,sol(nfit),serr(nfit,2),err(nfit),Psec,
     .  M1,M2,R1,aConst,asemi,Pi,Pid2,G,sb,Msun,Mearth,Mjup,Rjup,Rsun,
     .  Lsun,incl,R2,eccn,tPi,K,AU,b,rdrs,adrs,tdepth,time(nmax),
     .  flux(nmax),ferr(nmax),itime(nmax),MOSTtime,tdur,transitdur,
     .  time1(1),exptime1(1),dtype1(1),tmodel1(1),sn,epo,per,ph1,ph2,ph,
     .  pcut(nmax),std,stdev,pmean
      character*3 titles(nfit)
      character*80 inputsol,obsfile

      if(iargc().lt.2) goto 902
      
C     Constants
      Pi=acos(-1.d0)   !Pi
      Pid2=Pi/2.0d0
      tPi=Pi*2.0d0
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      sb=5.6704d-8 !W m^-2 K^-4 Stefan-Boltzman constant
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Mjup=317.833d0*Mearth !kg  mass of Jupiter  
      Rjup=142980.0d0*1000.0d0/2.0d0 !m  radius of Jupiter     
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Lsun=3.839d26 !W Solar Luminosity  
      AU=1.49598e11    
      
      call getarg(1,obsfile) !get filename for input data
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,flux,ferr,itime,MOSTtime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,MOSTtime)
      close(nunit)!release unit number as we are done with the file.
      if(npt.le.0)then
        write(6,*) -1.0
        goto 999
      endif
      
      call getarg(2,inputsol)
      nunit=10
      open(unit=nunit,file=inputsol,status='old',err=901)
      call getfitpars(nunit,nfitg,solg,serrg,Dpvaryg,errg,doe,toff)
      
      tdur=transitdur(nfit,solg)/86400.0d0 !transit duration in days
      
      npt1=1
      time1(1)=solg(7)
      exptime1(1)=1765.5/86400.0d0
      dtype1(1)=0
      call transitmodel(npt1,time1,exptime1,dtype1,tmodel1,nfit,solg)
      tdepth=(1.0d0-10**((tmodel1(1)-solg(8))/-2.5d0))
c      write(0,*) "TDEPTH:",tdepth
      
c      write(6,*) int(tdepth+0.5),tdur
      
      epo=solg(7)
      per=solg(5)
      ph1=(tdur)/solg(5)/2.0d0
      ph2=1.0d0+(-tdur)/solg(5)/2.0d0
c      write(6,*) ph1,ph2
      
      j=0
      do 10 i=1,npt
        ph=(time(i)-epo)/per-int((time(i)-epo)/per)
        if(ph.lt.0.0) ph=ph+1.0
        if((ph.gt.ph1).and.(ph.lt.ph2))then
            j=j+1
            pcut(j)=(1.0d0-10**((flux(i)-solg(8))/-2.5d0))
c            write(6,*) time(i),pcut(j)
        endif
 10   continue
 
      std=stdev(j,pcut,pmean)
      
      sn=tdepth/std*sqrt(dble(npt-j))
c      write(0,*) "STATS:", npt,std,j
      write(6,*) sn
 
c      
c      write(6,*) int(tdepth+0.5)
c 500  format(F6.4)
      
      
      goto 999
 901  write(0,*) "Cannot open" ,inputsol
      goto 999
 902  write(0,*) "Usage: transitsn photometry f1.dat"
      goto 999
 903  write(0,*) "Cannot open", obsfile
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function stdev(npt,pts,mean)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculates standard deviation of data set given the mean.
      implicit none

      integer npt,i
      double precision pts(npt),mean,s,ep,adev,p,var,sdev,mean2

      s=0.
      do 11 i=1,npt
         s=s+pts(i)
 11   continue
      mean=s/npt

      ep=0.
      var=0.
      do 10 i=1,npt
         s=pts(i)-mean
         ep=ep+s
         p=s*s
         var=var+p
 10   continue
      var=(var-ep**2/npt)/(npt-1)
      sdev=sqrt(var)

      stdev=sdev

      return
      end      