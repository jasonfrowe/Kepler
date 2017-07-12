      program transitplot3
      implicit none
      integer iargc,nunit,ntype,id,id2,kID,Teff,bins,nmax,nfit,npt,i,
     .  nsensor
      parameter(nmax=650000,nfit=18)
      integer dtype1(nmax),oddeven(nmax),flag(nmax)
      real px(nmax),py(nmax),pz(nmax)
      double precision Kmag,logg,rad,eoff,escale,xwidth,time(nmax),
     .  mag(nmax),merr(nmax),itime(nmax),Zerotime,sol(nfit),
     .  serr(nfit,2),Dpvary(nfit),err(nfit),doe,toff,zpt,flux(nmax),
     .  ferr(nmax),bb(4),phase(nmax),bbp(4),bbp2(4),avgitime,
     .  avgvtime,mean,bbp3(4),transitdur,tdur,ph1,ph2,sym,tmodel(nmax),
     .  tdepth,dtype(1),timeone(1),exptime(1),flux2(nmax),ferr2(nmax),
     .  kerr,phase2(nmax),pts(nmax),Thpars(3,2)
      character*80 parsfile,filename,fitparsfile

      nsensor=0 !sensor text if ==1, no text if ==2

      if(iargc().lt.1) goto 901 !number of required commandline pars
      call getarg(1,parsfile) !get first parameter - filename

      nunit=10
      open(unit=nunit,file=parsfile,status='old',err=903)
      call getplotpars(nunit,filename,ntype,fitparsfile,id,id2,kID,
     .  Kmag,Teff,logg,rad,eoff,escale,xwidth,bins)
      close(nunit)
      
      write(6,*) "KID",kID
      
      open(unit=nunit,file=filename,status='old',err=902)
      
      if (ntype.eq.0)then
        call readdata(nunit,nmax,npt,time,mag,merr,itime,Zerotime)
      elseif(ntype.eq.1)then
        Zerotime=54900.0
        call readkeplc(nunit,nmax,npt,time,mag,merr,itime,Zerotime)
      endif
      close(nunit)

      if(kmag.gt.0.0)then !estimate of expected scatter
C       quadratic version
c        kerr=(4.0d0*(kmag-8.0d0)**2.0+10.0d0)*1.0d-6
C       Powerlaw version
         kerr=(10.0d0**(0.23*kmag-0.9))*1.0d-6
      else
        kerr=100.0*1.0d-6
      endif 
      kerr=2.5*log10(1.0d0+kerr) !convert flux error to mag
      
      do 11 i=1,npt
        merr(i)=kerr
        dtype1(i)=0
 11   continue
 
      open(unit=nunit,file=fitparsfile,status='old',err=904)
      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff,eoff)
      close(nunit)

      
c      if(bins.gt.0) bins=int(sol(5)*1440.0d0/30.0d0)
c      write(6,*) "bins",bins

      zpt=sol(8)

      do 10 i=1,npt
        dtype1(i)=0
        mag(i)=mag(i)-zpt
 10   continue
      sol(8)=0.0
      
C     Convert magnitudes to flux
      call mag2flux(npt,mag,flux)
      call mag2flux(npt,merr,ferr)
      do 14 i=1,npt
        ferr(i)=1.0d0-ferr(i)
c        write(6,*) time(i),flux(i),ferr(i)
 14   continue
 

      toff=0.75-(sol(7)/sol(5)-int(sol(7)/sol(5)))
      if(toff.lt.0.0)toff=toff+1.0
      tdur=transitdur(nfit,sol)/86400.0d0
      
      call transitstats(npt,time,phase,flux,ferr,itime,dtype1,pts,
     .  tmodel,nfit,sol,toff,Teff,err,Thpars)
      
      timeone(1)=sol(7)
      dtype(1)=0
      exptime(1)=1765.5/86400.0d0
      call transitmodel(1,timeone,exptime,dtype,tmodel,nfit,sol)
      tdepth=(1.0d0-10**((tmodel(1)-sol(8))/-2.5d0))
      write(6,*) "Tdepth: ",tdepth*1.0d6
      
c      tdur=0.30
      write(0,*) tdur,sol(5),tdur/sol(5)
      ph1=0.75-1.5d0*tdur/sol(5)
      if(ph1.lt.0.5)ph1=0.5
      ph2=0.75+1.5d0*tdur/sol(5)
      if(ph2.gt.1.0)ph2=1.0
 
      call opengraphics()
      call pgslw(2)
      
      call plotflux(npt,time,flux,ferr,-1.0d0,-1.0d0,1,Zerotime,nmax,
     .  px,py,pz,bb,nfit,sol,id,id2,nsensor,kID,Kmag,Teff,logg,rad,
     .  tdepth)
      bins=100
      call plotphflux(npt,time,flux,ferr,sol(5),0.0d0,1.0d0,2,sym,nmax,
     .  px,py,pz,phase,toff,1,bbp,tdepth,flux2,ferr2,bins,nfit,sol,err,
     .  nsensor,thpars,flag,rad,tmodel,itime,dtype1,tdur)
      bins=int(sol(5)*1440.0d0/30.0d0)
      call plotphflux(npt,time,flux,ferr,sol(5),ph1,ph2,3,sym,nmax,
     .  px,py,pz,phase,toff,1,bbp2,tdepth,flux2,ferr2,bins,nfit,sol,err,
     .  nsensor,thpars,flag,rad,tmodel,itime,dtype1,tdur)
      call plotoddeven(npt,time,flux,ferr,sol(5),ph1,ph2,4,sym,nmax,
     .  px,py,pz,phase,toff,1,bbp3,tdepth,flux2,ferr2,bins,oddeven,
     .  phase2,doe,nsensor,flag)
      avgitime=mean(npt,itime)
c      avgvtime=mean(nptv,vetime)
      call modelplotflux(nfit,sol,nmax,phase,bb,bbp,bbp2,bbp3,
     .  avgitime,avgvtime,toff)
      
      call pgslw(1)
      call closegraphics()

      goto 999
 901  write(6,*) "Usage: transitplot <parsfile>"
      goto 999
 902  write(6,*) "Cannot open ",filename
      goto 999
 903  write(6,*) "Cannot open ",parsfile
      goto 999
 904  write(6,*) "Cannot open ",fitparsfile
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
      subroutine transitstats(npt,time,phase,flux,ferr,itime,dtype,pts,
     .  tmodel,nfit,sol,toff,Teff,err,Thpars)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nfit,i,j,Teff,dtype(npt)
      double precision time(npt),flux(npt),tmodel(npt),sol(nfit),
     .  itime(npt),lchi,tflux,mchi,zpt,pdur,transitdur,Psec,Pi,tPi,
     .  tcentre,ecentre,toff,phase(npt),pmin,pmax,G,Msun,Mearth,
     .  Mjup,M1,M2,asemi,AU,Eanom,eccn,Manom,phi,trueanomaly,w,Ab,Rjup,
     .  mean1,mean2,std1,std2,pts(npt),pmin2,pmax2,stdev,Rsun,
     .  Teq,Ag,f,sb,Lstar,Lp,Tp,Lsun,Agerr,err(nfit),asemierr,aConst,
     .  Teqerr,Tefferr,R1,Lstarerr,Lperr,temp(4),Tperr,Thpars(3,2),
     .  R1err,ferr(npt),tidal,K,incl,Pid2,M1err

      Pi=acos(-1.d0)   !Pi
      Pid2=Pi/2.0d0
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      sb=5.6704d-8 !W m^-2 K^-4 Stefan-Boltzman constant
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Mjup=317.833d0*Mearth !kg  mass of Jupiter  
      Rjup=142980.0d0*1000.0d0/2.0d0 !m  radius of Jupiter     
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Lsun=3.839d26 !W Solar Luminosity
      Psec=sol(5)*8.64d4 !sec ; period of planet
      M1=sol(1)*Msun !kg ; mass of star
      M2=sol(2)*Mjup !kg ; mass of planet
      R1=sol(3)*Rsun !m ; stellar radius
      AU=1.49598e11
      incl=sol(6)
      if(incl.gt.90.0d0)incl=180.0-incl     
      incl=Pi*(90.0d0-incl)/180.0d0
      
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      asemi=(M1+M2)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst  
C     semi-major axis error needs to include mass error!!
c      asemierr=2.0d0*aConst*Psec**(-1.0d0/3.0d0)*err(5)*8.64d4/3.0d0
      M1err=0.05*M1!err(1)*Msun!0.3*M1
      temp(1)=M1**(-2.0d0/3.0d0)*(M1err)/3.0 !assume 30% error in mass
      temp(2)=2.0d0*Psec**(-1.0d0/3.0d0)*err(5)*8.64d4/3.0d0
      asemierr=asemi*asemi*( (temp(1)/M1**(1.0d0/3.0d0))**2.0d0 + 
     .  (temp(2)/Psec**(2.0d0/3.0d0))**2.0d0  )
      asemierr=sqrt(asemierr)
      write(6,501) "Asemi ",asemi/AU,"+/-",asemierr/AU
 501  format(A6,F8.4,A3,F8.4)
C     Ellipsoidal amplitude
      tidal=M2/M1*(R1/asemi)**3.0d0*1.0d6
      write(6,505) "Ellip_amp: ",tidal
 505  format(A11,F7.4)
C     Calculated measured albedo
      Ag=2.0d0*sol(16)*1.0d-6*asemi*asemi/
     .  (3.0d0*sol(4)*sol(4)*Rjup*Rjup)
      Agerr=Ag*Ag*( (err(16)/sol(16))**2.0d0 + 
     .  (2.0d0*asemierr/asemi)**2.0d0 +
     .  (2.0d0*err(4)/sol(4))**2.0d0 )
      Agerr=sqrt(Agerr)
      write(6,500) "Ag: ",Ag,"+/-",Agerr
      write(6,504) "AG_Fp/F*: ",
     .  Ag*sol(4)*sol(4)*Rjup*Rjup/(asemi*asemi)*1.0d6
 504  format(A10,F7.3)
 500  format(A4,F7.3,A3,F7.3)
      Ab=0.3 !assumed bond albedo for planet temperature
      f=1.0 !circulation factor
C     Calculate equilibrium temperature based on assumed albedo and f
      Teq=dble(Teff)*sqrt(sol(3)*Rsun/(2.0d0*asemi))*
     .  (f*(1-Ab))**(0.25d0)
      Tefferr=150.0 !Assumed error in the stellar effective temperature
      R1err=R1*0.3!err(3)*Rsun!R1*0.3 !let assume a 30% error in stellar radius
      Teqerr=Teq*Teq*(
     .  (Tefferr/dble(Teff))**2.0d0 +
     .  ((R1err/R1)**2.0d0)/4.0d0 +
     .  ((asemierr/asemi)**2.0d0)/4.0d0)
      Teqerr=sqrt(Teqerr)
      write(6,502) "Teq: ",int(Teq+0.5d0),"+/-",int(Teqerr+0.5d0)
 502  format(A5,I5,A3,I5)
C     Calculate luminosity of the Star
      Lstar=4.0*Pi*(sol(3)*sol(3)*Rsun*Rsun)*sb*(dble(Teff))**4.0
      Lstarerr=Lstar*Lstar*( (R1err/R1)**2.0d0 + 
     .  (4.0d0*Tefferr/dble(Teff))**2.0d0 )  !assume 30% error in radius
      Lstarerr=sqrt(Lstarerr)
c      Lstarerr=2.4*Lsun
C     Calculate the luminosity of the planet
      Lp=Lstar*sol(16)*1.0d-6!*2.0d0
      Lperr=Lp*Lp*( (Lstarerr/Lstar)**2.0 + (err(16)/sol(16))**2.0 )
      Lperr=sqrt(Lperr)
      write(6,503) "Ls: ",Lstar/Lsun,"+/-",Lstarerr/Lsun
      write(6,503) "Lp: ",Lp/Lsun,"+/-",Lperr/Lsun
 503  format(A4,F11.5,A3,F11.5)
C     Measured Temperature of the planet
      if(Lp.le.0.0)then
        Tp=0.0
        Tperr=0.0
      else
        temp(1)=Lp/(4.0*Pi*Rjup*Rjup*sol(4)*sol(4)*sb)
        Tp=(temp(1))**0.25
        temp(2)=(temp(1)**(-3.0d0/4.0d0))/4.0d0
        temp(3)=2.0d0*err(4)/sol(4)
        temp(4)=temp(1)*temp(1)*( (Lperr/Lp)**2.0 + 
     .    (temp(3)*temp(3)/(4.0d0*Pi*Rjup*Rjup*sol(4)*sol(4)*sb))**2.0 )
        temp(4)=sqrt(temp(4))
        Tperr=temp(2)*temp(4)
      endif
      
      write(6,502) "Tp:  ",int(Tp+0.5),"+/-",int(Tperr+0.5)
      
C     Save thermal parameters
      Thpars(1,1)=Ag
      Thpars(1,2)=Agerr
      Thpars(2,1)=Teq
      Thpars(2,2)=Teqerr
      Thpars(3,1)=Tp
      Thpars(3,2)=Tperr
      
c      eccn=sol(14) !the eccentricity of the planet orbit
      eccn=sqrt(sol(14)*sol(14)+sol(15)*sol(15)) !eccentricity

C     Estimate of radial velocity amplitude      
      K=2.0*pi*G*M2**3*(sin(incl+Pid2))**3/
     .  (Psec*(1.0d0-eccn*eccn)**(3.0d0/2.0d0)*(M1+M2)*(M1+M2))
      K=K**(1.0d0/3.0d0)
      write(6,*) "K: ",K
      
      write(6,*) "Doppler: ",6.5*K/2.99792458e8
      write(6,*) "Tides : ",3.37*M2/M1*(R1/asemi)**3*(sin(incl+Pid2))**3
      
      call transitmodel(npt,time,itime,dtype,tmodel,nfit,sol)
      call phasept(npt,time,phase,sol(5),toff)
      
C     Get centre of transit and eclipse 
      Pi=acos(-1.d0)   
      tPi=2.0d0*Pi  
      
C     Find the position of the transit and eclipse
      if(eccn.eq.0.0)then 
        phi=pi-2.0*pi*(sol(7)/sol(5)-int(sol(7)/sol(5)))
        if(phi.lt.0.0) phi=phi+tpi
C       The case of the circular orbit it trivial
        tcentre=0.5-phi/tpi+toff !find center of transit phase
        if(tcentre.lt.0.0d0) tcentre=tcentre+1.0d0!make sure between 0-1
        if(tcentre.gt.1.0d0) tcentre=tcentre-1.0d0
        ecentre=1.0-phi/tpi+toff !find center of eclipse phase
        if(ecentre.lt.0.0d0) ecentre=ecentre+1.0d0!make sure between 0-1
        if(ecentre.gt.1.0d0) ecentre=ecentre-1.0d0
      else
c        w=Pi*sol(15)/180.0d0!+Pi
        w=acos(sol(14)/eccn)
        Eanom=trueanomaly(eccn,w)
        Manom=Eanom
        call invKepler(Eanom,Manom,eccn)
        phi=sol(7)+Manom!-Pi
        write(0,*) sol(7),Manom,w

        tcentre=0.5-phi/tpi+toff
        ecentre=phi/tpi+toff
c        tcentre=0.5-phi/tpi+toff !find center of transit phase
c        ecentre=1.0-phi/tpi+toff !find center of eclipse phase
        tcentre=tcentre-int(tcentre)
        if(tcentre.lt.0.0d0)tcentre=tcentre+1.0
        ecentre=ecentre-int(ecentre)
        if(ecentre.lt.0.0d0)ecentre=ecentre+1.0
        
        write(0,*) "t:",tcentre,ecentre

c        arad(j)=distance(asemi,eccn,Tanom(j))
c        x2(j)=arad(j)*Sin(Tanom(j)+w)
c        y2(j)=arad(j)*Sin(incl)*Cos(Tanom(j)+w)
      endif
      
      Psec=sol(5)*24.0*60.0*60.0 !sec ; period of planet
      pdur=transitdur(nfit,sol)/Psec/2.0d0 !transit duration
      
C     calculate chi-squared with straight line (lchi) and model (mchi)
      lchi=0.0d0
      mchi=0.0
      zpt=real(10.0**(sol(8)/-2.5))
      pmin=tcentre-pdur
      pmax=tcentre+pdur
      do 10 i=1,npt
        if((phase(i).ge.pmin).and.(phase(i).le.pmax))then
            tflux=10.0**(tmodel(i)/-2.5)
            mchi=mchi+(flux(i)-tflux)*(flux(i)-tflux)/(ferr(i)*ferr(i))
            lchi=lchi+(flux(i)-zpt)*(flux(i)-zpt)/(ferr(i)*ferr(i))
        endif
 10   continue
c      mchi=mchi/dble(npt-1)
c      lchi=lchi/dble(npt-1)
 
      write(6,*) "Transit:",lchi-mchi
      
C     Do the same for the eclipse
      lchi=0.0d0
      mchi=0.0d0
      zpt=real(10.0**(sol(8)/-2.5))
      pmin=ecentre-2.0d0*pdur
      pmax=ecentre+2.0d0*pdur
      j=0
      do 11 i=1,npt
        if((phase(i).ge.pmin).and.(phase(i).le.pmax))then
            tflux=10.0**(tmodel(i)/-2.5)
            mchi=mchi+(flux(i)-tflux)*(flux(i)-tflux)/(ferr(i)*ferr(i))
            lchi=lchi+(flux(i)-zpt)*(flux(i)-zpt)/(ferr(i)*ferr(i))
            j=j+1
c            write(6,*) j,zpt,flux(i),tflux
c            read(5,*)
        endif
 11   continue
  
      write(6,*) "Eclipse:",lchi-mchi
      
      pmin=ecentre-pdur
      pmax=ecentre+pdur
      j=0
      do 12 i=1,npt
        if((phase(i).ge.pmin).and.(phase(i).le.pmax))then
            j=j+1
            pts(j)=flux(i)
        endif
 12   continue
      std1=stdev(j,pts,mean1)
      
      pmin2=ecentre-2.0d0*pdur
      pmax2=ecentre+2.0d0*pdur
      j=0
      do 13 i=1,npt
        if(((phase(i).ge.pmin2).and.(phase(i).le.pmin)).or.
     .    ((phase(i).ge.pmax).and.(phase(i).le.pmax2)))then
            j=j+1
            pts(j)=flux(i)
        endif
 13   continue
      std2=stdev(j,pts,mean2)
      
      write(6,*) mean1,std1
      write(6,*) mean2,std2 
      
             
      return
      end
