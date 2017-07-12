CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program transitplot
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Plots for TCERT
      implicit none
      integer iargc,nunit,nmax,npt,ntype,nfit,id,kID,Teff,i,id2,bins
      parameter(nmax=1200000,nfit=18)
      integer sym(nmax),ps(nmax),nsensor,dtype1(nmax),flag(nmax)
      real px(nmax),py(nmax),pz(nmax)
      double precision time(nmax),mag(nmax),merr(nmax),itime(nmax),
     .  Zerotime,sol(nfit),serr(nfit,2),Dpvary(nfit),toff,flux(nmax),
     .  Kmag,logg,rad,eoff,phase(nmax),escale,xwidth,transitdur,
     .  tmodel(nmax),pts(nmax),zpt,ferr(nmax),kerr,err(nfit),doe,
     .  Thpars(3,2),phase2(nmax),flux2(nmax),ferr2(nmax)
C     RV-related varibles
      integer nptv
      double precision vtime(nmax),vel(nmax),verr(nmax),vetime(nmax)
      character*80 filename,parsfile,fitparsfile
      
      nsensor=0!sensor text if ==1, no text if ==2
      bins=0 !if you want phase binning (use something like 200)
            
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

c      if(kmag.gt.0.0)then !estimate of expected scatter
cC       quadratic version
cc        kerr=(4.0d0*(kmag-8.0d0)**2.0+10.0d0)*1.0d-6
cC       Powerlaw version
c         kerr=(10.0d0**(0.23*kmag-0.9))*1.0d-6
c      else
c        kerr=100.0*1.0d-6
c      endif 
c      kerr=2.5*log10(1.0d0+kerr) !convert flux error to mag
c      write(0,*) "KMAG,KERR:",kmag,kerr
      do 11 i=1,npt
c        merr(i)=kerr
        dtype1(i)=0
c        itime(i)=10.0/(24.0*60.0)
 11   continue
            
      open(unit=nunit,file=fitparsfile,status='old',err=904)
      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff,eoff)
      close(nunit)
      
      if(bins.gt.0) bins=int(sol(5)*1440.0d0/30.0d0)
      write(6,*) "bins",bins
      
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

cC     Now lets get the RV data
c      if(rvfile.eq.'null')then
c        nptv=0
c      else
c        nunit=10
c        open(unit=nunit,file=rvfile,status='old',err=905)
c        call readrv(nunit,nmax,nptv,vtime,vel,verr,vetime)
c        close(nunit)
cc        write(6,*) "rv:",(vtime(i),i=1,nptv)
c      endif
 
      call pgopen('?')
      call pgpap(8.0,0.85)
c      call pgask(.false.)
c      call pgsch(0.8)
      call pgsubp(1,2)
      
      write(6,*) "Transit duration (h)",transitdur(nfit,sol)/3600.0
      if(xwidth.le.0.0) xwidth=3.0*transitdur(nfit,sol)/3600.0
      call transitstats(npt,time,phase,flux,ferr,itime,dtype1,pts,
     .  tmodel,nfit,sol,toff,Teff,err,Thpars)
      
      call pgvport(0.15,0.85,0.2,0.9)
      call panelA(npt,px,py,pz,time,flux,ferr,1,1,id,id2,kID,Kmag,Teff,
     .  logg,rad,Zerotime,nsensor,nfit,sol)
      call pgvport(0.15,0.85,0.2,1.0)
      call pgvport(0.15,0.85,0.2,0.9)
      call panelB(nmax,npt,px,py,pz,time,flux,ferr,itime,phase,1,2,nfit,
     .  sol,toff,eoff,escale,sym,ps,xwidth,Zerotime,err,doe,Thpars,
     .  nsensor,bins,phase2,flux2,ferr2,flag)

      call pgclos()
      
      goto 999
 901  write(6,*) "Usage: transitplot <parsfile>"
      goto 999
 902  write(6,*) "Cannot open ",filename
      goto 999
 903  write(6,*) "Cannot open ",parsfile
      goto 999
 904  write(6,*) "Cannot open ",fitparsfile
      goto 999
c 905  write(6,*) "Cannot open ",rvfile
c      goto 999
 999  end

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
      subroutine panelC(nmax,npt,px,py,pz,time,flux,ferr,itime,phase,
     .  npx,npy,nfit,sol,toff,eoff,escale,sym,ps,xwidth,Zerotime,err,
     .  doe,Thpars,nsensor,bins,phase2,flux2,ferr2,flag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,npx,npy,nfit,nplot,i,nmax,np,sym(nmax),ps(nmax),bins,
     .  npt2,nplot2,nsensor
      parameter(nplot=100000)
      integer dtype(nplot),flag(npt)
      real px(nmax),py(nmax),rbb(4),rescale,zpt,bscale,offset,pz(nmax)
      double precision time(npt),flux(npt),sol(nfit),ptime(nplot),
     .  itime(npt),dnpt,bb(4),tmodel(nplot),eoff,phase(nmax),toff,
     .  escale,xwidth,tcentre,ecentre,Pi,tPi,pitime(nplot),itmean,
     .  epoch,Zerotime,mintime,poff,Eanom,Manom,Tanom,eccn,w,
     .  trueanomaly,phi(2),ferr(nmax),temp,err(nfit),redj,eoff2,doe,
     .  Thpars(3,2),terr(nplot),phase2(nmax),flux2(nmax),ferr2(nmax)
      character*80 label

      if(xwidth.gt.1.0d0) xwidth=min(xwidth/(sol(5)*24.0),0.25)

      redj=0.09205 !radius of Earth divided by Jupiter

      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi


      eccn=sqrt(sol(14)*sol(14)+sol(15)*sol(15)) !eccentricity

C       The case of the circular orbit is trivial
        tcentre=sol(7)/sol(5)-floor(sol(7)/sol(5))
        if(tcentre.lt.0.0d0) tcentre=tcentre+1.0d0!make sure between 0-1
        if(tcentre.gt.1.0d0) tcentre=tcentre-1.0d0
        ecentre=sol(7)/sol(5)-floor(sol(7)/sol(5))+0.5
        if(ecentre.lt.0.0d0) ecentre=ecentre+1.0d0!make sure between 0-1
        if(ecentre.gt.1.0d0) ecentre=ecentre-1.0d0

c      write(0,*) "tc,ec:",tcentre,ecentre

C     convert escale to real (used for eclipse scale)
      rescale=real(escale)

      call phasept(npt,time,phase,sol(5),0.0d0)
      do 5 i=1,npt
        phase(i)=phase(i)-tcentre
        if(phase(i).lt.-0.5d0) phase(i)=phase(i)+1.0d0
        if(phase(i).gt. 0.5d0) phase(i)=phase(i)-1.0d0
 5    continue

c      if(bins.gt.0) call phaseptsp(npt,time,phase,sol(5),0.0d0,sym)
      do 8 i=1,npt !copy data, because binning is destructive
        phase2(i)=phase(i)
        flux2(i)=flux(i)
        ferr2(i)=ferr(i)
 8    continue
      npt2=npt
c      write(0,*) "hello1"
      if(bins.gt.0) call binp(npt2,phase2,flux2,ferr2,bins,flag)!binning

c      write(0,*) "hello2"

C     Find means of dataset
      itmean=0.0d0 !itmean is the mean exposure time
      bb(1)=time(1) !lower x-axis bound
      bb(2)=time(1) !high x-axis bound
      do 11 i=1,npt
        bb(1)=min(time(i),bb(1)) !mix min and max times
        bb(2)=max(time(i),bb(2))
        itmean=itmean+itime(i)
 11   continue
      itmean=itmean/dble(npt) !mean exposure time
C     time of first observation - used to get the epoch below
      mintime=bb(1)

C     Find bounds of data that will be plotted.
      bb(3)= 99.9d30!flux(1) !lower y-axis bound
      bb(4)=-99.9d30!flux(1) !high y-axis bound
      eoff2=-99.9d30!find max during secondary
      do 9 i=1,npt2
        if(abs(phase2(i)).lt.xwidth)then
            bb(3)=min(flux2(i),bb(3))
            bb(4)=max(flux2(i),bb(4))
        endif
        if((phase2(i).gt.0.5-xwidth).or.(phase2(i).lt.-0.5+xwidth))then
            eoff2=max(eoff2,flux2(i)-1.0)
        endif
 9    continue
      if(eoff.le.0.0) eoff=2.0*abs(1.0-bb(4)) !auto offsets
c      eoff=2.0*abs(1.0-bb(4))  !remove this line
c      eoff=2.0*abs(1.0-bb(4))
      if(eoff2.le.0.0) eoff2=0.0
C     Choose between next two line for top of plot blank space
c      eoff2=eoff2+eoff2*0.15 !tag on 15% for nice plot boundaries
      eoff2=eoff
      if(nsensor.eq.2) then !when labels go away, more space for play
        bb(3)=bb(3)-(1.0-bb(3))*0.05
      else
        bb(3)=bb(3)-(1.0-bb(3))*0.25
      endif

      dnpt=dble(nplot-1) !convert integer to real for division
      do 10 i=1,nplot  !generate 1 phase of model data
C     Psec contains the orbital Period of the planet in seconds
        ptime(i)=bb(1)+(bb(2)-bb(1))*dble(i-1)/dnpt !observation times
        if(bins.gt.0)then
            pitime(i)=sol(5)/dble(bins)
        else
            pitime(i)=itmean !mean integration times
        endif
        dtype(i)=0
 10   continue

      call pgpanl(npx,npy)
      bb(1)=-real(xwidth)
      bb(2)= real(xwidth)
      rbb(1)=real(bb(1))
      rbb(2)=real(bb(2))
      rbb(3)=real(bb(3)-0.35*(bb(4)-bb(3)))
      rbb(4)=real(eoff+eoff2)+1.0

      call pgwindow(24.0*rbb(1)*real(sol(5)),24.0*rbb(2)*real(sol(5)),
     .  rbb(3),rbb(4))

      call pgsch(1.5) !2.0
      if(nsensor.le.2) call pgslw(3)
cc      call pgbox('BCNTS1',0.0,0,'BNTSV1',0.0,0)
c      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgsch(2.0) !c1 made bigger labels.
      write(label,503) "Phase (Hours)"
 503  format(A13)
      call pgptxt((rbb(1)+rbb(2))/2.0,
     .  rbb(3)-0.14*(rbb(4)-rbb(3)),
     .  0.0,0.5,label)
cc      call pgptxt(rbb(1)-0.11*(rbb(2)-rbb(1)),(rbb(4)+rbb(3))/2.0,
cc     .  90.0,0.5,"Relative Flux")
      call pgslw(1)
      call pgsch(1.0)

      call phaseptsp(npt,time,phase,sol(5),0.0d0,sym)
      if(bins.gt.0) call binp(npt,phase,flux,ferr,bins,flag)

      np=0
      do 14 i=1,npt
        temp=phase(i)-tcentre
        if(temp.lt.-0.5) temp=temp+1.0
        if(temp.gt. 0.5) temp=temp-1.0
        if((temp.ge.bb(1)).and.
     .    (temp.le.bb(2)))then
            np=np+1
            px(np)=real(temp)!real(phase(i))-tcentre
            py(np)=real(flux(i))
            pz(np)=real(ferr(i))
            ps(np)=sym(i)
c            write(6,*) abs(px(np)),py(np),pz(np)
        endif
 14   continue
      call pgsch(2.0)
c      call pgerrb(6,np,px,py,pz,1.0)
cc      call pgpnts(np,px,py,ps,np)
      call pgsch(1.0)

      zpt=real(10.0**(sol(8)/-2.5))
      bscale=(rbb(4)-rbb(3))/(2.0*rescale)
      offset=((zpt+eoff)-(rbb(4)+rbb(3))/2.0)/rescale
      call pgwindow(rbb(1),rbb(2),zpt-bscale,zpt+bscale)
      call pgsch(2.0) !2.0
      if(nsensor.le.2) call pgslw(3)
      call pgptxt(rbb(1)-0.12*(rbb(2)-rbb(1)),zpt,
     .  90.0,0.5,"Folded Amplitude (ppm)")
      call pgslw(1)
      call pgsch(1.0)
c      call pgbox('',0.0,0,'CMTSV1',0.0,0)
      np=0
      do 15 i=1,npt
        temp=phase(i)-ecentre
        if(temp.lt.-0.5) temp=temp+1.0
        if(temp.gt. 0.5) temp=temp-1.0
        if((temp.ge.bb(1)).and.
     .    temp.le.bb(2))then
            np=np+1
            px(np)=real(temp)!real(phase(i))-ecentre
            py(np)=real(flux(i))+offset
            pz(np)=real(ferr(i))
        endif
 15   continue
c      call pgerrb(6,np,px,py,pz,1.0)
      call pgpt(np,px,py,4)
c      write(6,*) 1.0-(zpt-bscale-offset),1.0-(zpt+bscale-offset)
      call pgwindow(rbb(1),rbb(2),1.0e6*(1.0-(zpt-bscale-offset)),
     .  1.0e6*(1.0-(zpt+bscale-offset)))
      call pgsch(1.5)!2.0
      if(nsensor.le.2) call pgslw(3)
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
cc      call pgbox('',0.0,0,'CMTSV1',0.0,0)
c1      call pgbox('',0.0,0,'CTSV1',0.0,0)
      call pgsch(1.0)
      call pgslw(1)
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))

      call transitmodel(nplot,ptime,pitime,dtype,tmodel,nfit,sol)
      call phasept(nplot,ptime,phase,sol(5),0.0d0)
      do 16 i=1,nplot
        terr(i)=1.0
 16   continue

      nplot2=nplot
c      call binp(nplot2,phase,tmodel,terr,bins,flag)

      do 6 i=1,nplot2
        phase(i)=phase(i)-tcentre
        if(phase(i).lt.-0.5) phase(i)=phase(i)+1.0
        if(phase(i).gt. 0.5) phase(i)=phase(i)-1.0
 6    continue
      call sort2(nplot2,phase,tmodel) !sort in phase

      np=0
      do 12 i=1,nplot2
        temp=phase(i)
        if((temp.ge.bb(1)).and.
     .    (temp.le.bb(2)))then
            np=np+1
            px(np)=real(temp)!real(phase(i))-tcentre
            py(np)=real(10.0**(tmodel(i)/-2.5))
        endif
 12   continue

      call pgsch(2.0)
      call pgslw(3) !thicker lines
      call pgsci(2)
c      call pgline(np,px,py)

      call pgwindow(rbb(1),rbb(2),zpt-bscale,zpt+bscale)

      do 7 i=1,nplot2
        phase(i)=phase(i)+tcentre-ecentre
        if(phase(i).lt.-0.5) phase(i)=phase(i)+1.0
        if(phase(i).gt. 0.5) phase(i)=phase(i)-1.0
 7    continue
      call sort2(nplot2,phase,tmodel) !sort in phase

      np=0
      do 13 i=1,nplot2
        if((phase(i).ge.bb(1)).and.
     .    (phase(i).le.bb(2)))then
            np=np+1
            px(np)=real(phase(i))
            py(np)=real(10.0**(tmodel(i)/-2.5))+offset
        endif
 13   continue
      call pgsci(3)
      call pgline(np,px,py)

      call pgsci(1.0)
      call pgslw(1)
      call pgsch(1.0)

C     Add labels.
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgsch(1.5)
      if(nsensor.le.2) call pgslw(3) !thicker lines
      write(label,500) "incl=",sol(6)
 500  format(A5,F5.2)
      write(6,*) "NSENSOR",nsensor
      if(nsensor.lt.2) call pgptxt(rbb(2)-0.03*(rbb(2)-rbb(1)),
     .  rbb(3)+0.055*(rbb(4)-rbb(3)),
     .  0.0,1.0,label)
      if(sol(4)/redj.lt.3.0d0)then
        write(label,501) "Rp=",sol(4)/redj,"Re"
      else
        write(label,501) "Rp=",sol(4),"Rj"
      endif
 501  format(A3,F5.2,1X,A2)
      if(nsensor.lt.2) call pgptxt(rbb(2)-0.03*(rbb(2)-rbb(1)),
     .  rbb(3)+0.11*(rbb(4)-rbb(3)),
     .  0.0,1.0,label)
      if(err(16).gt.0.0d0)then
        if(sol(16)/err(16).lt.10.0)then
            write(label,505) "Eclip-sig ",sol(16)/err(16)
        else
            write(label,512) "Eclip-sig ",sol(16)/err(16)
        endif
 505    format(A10,F4.1)
 512    format(A10,F5.1)
        if(nsensor.lt.2) call pgptxt(rbb(2)-0.03*(rbb(2)-rbb(1)),
     .      rbb(3)+0.275*(rbb(4)-rbb(3)),
     .      0.0,1.0,label)
      endif
      write(6,*) "doe:",doe
      if(doe.gt.0.0d0)then
        write(label,505) "Depth-sig ",doe
        if(nsensor.lt.2) call pgptxt(rbb(2)-0.03*(rbb(2)-rbb(1)),
     .      rbb(3)+0.220*(rbb(4)-rbb(3)),
     .      0.0,1.0,label)
      endif


      poff=0.5-sol(7)/tpi

      epoch=sol(7)
      write(label,502) "Tc(E)=",epoch,"(BJD)+E(",sol(5),"days)"
 502  format(A6,F8.4,A8,F9.5,A5)
      if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .  rbb(3)+0.11*(rbb(4)-rbb(3)),
     .  0.0,0.0,label)
      write(label,504) "       ",err(7),"         ",err(5),"     "
 504  format(A7,F8.4,A9,F9.5,A5)
      if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .  rbb(3)+0.055*(rbb(4)-rbb(3)),
     .  0.0,0.0,label)

      if(sol(16)/err(16).ge.2.0d0)then !only plot if a reasonable detect
        if(Thpars(1,1).lt.10.0)then
            write(label,506) "Ag = ",Thpars(1,1),"(",Thpars(1,2),")"
 506        format(A5,F4.2,A1,F4.2,A1)
        else
            write(label,509) "Ag =",Thpars(1,1),"(",Thpars(1,2),")"
 509        format(A4,F6.2,A1,F6.2,A1)
        endif
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .      rbb(3)+0.275*(rbb(4)-rbb(3)),
     .      0.0,0.0,label)

        if(Thpars(2,1).lt.1.0d4)then
            write(label,507) "Teq=",int(Thpars(2,1)+0.5),"(",
     .          int(Thpars(2,2)+0.5),")"
 507        format(A4,I4,A1,I4,A1)
        else
            write(label,510) "Teq=",int(Thpars(2,1)+0.5),"(",
     .          int(Thpars(2,2)+0.5),")"
 510        format(A4,I5,A1,I5,A1)
        endif
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .      rbb(3)+0.220*(rbb(4)-rbb(3)),
     .      0.0,0.0,label)
        if(Thpars(3,1).lt.1.0d4)then
            write(label,508) "Teff=",int(Thpars(3,1)+0.5),"(",
     .          int(Thpars(3,2)+0.5),")"
 508        format(A5,I4,A1,I4,A1)
        else
            write(label,511) "Teff=",int(Thpars(3,1)+0.5),"(",
     .          int(Thpars(3,2)+0.5),")"
 511        format(A5,I5,A1,I5,A1)
        endif
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .      rbb(3)+0.165*(rbb(4)-rbb(3)),
     .      0.0,0.0,label)
      endif

      call pgsch(1.0) !return character and line widths to normal
      call pgslw(1)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine panelB(nmax,npt,px,py,pz,time,flux,ferr,itime,phase,
     .  npx,npy,nfit,sol,toff,eoff,escale,sym,ps,xwidth,Zerotime,err,
     .  doe,Thpars,nsensor,bins,phase2,flux2,ferr2,flag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,npx,npy,nfit,nplot,i,nmax,np,sym(nmax),ps(nmax),bins,
     .  npt2,nplot2,nsensor
      parameter(nplot=100000)
      integer dtype(nplot),flag(npt)
      real px(nmax),py(nmax),rbb(4),rescale,zpt,bscale,offset,pz(nmax)
      double precision time(npt),flux(npt),sol(nfit),ptime(nplot),
     .  itime(npt),dnpt,bb(4),tmodel(nplot),eoff,phase(nmax),toff,
     .  escale,xwidth,tcentre,ecentre,Pi,tPi,pitime(nplot),itmean,
     .  epoch,Zerotime,mintime,poff,Eanom,Manom,Tanom,eccn,w,
     .  trueanomaly,phi(2),ferr(nmax),temp,err(nfit),redj,eoff2,doe,
     .  Thpars(3,2),terr(nplot),phase2(nmax),flux2(nmax),ferr2(nmax)
      character*80 label
      
      if(xwidth.gt.1.0d0) xwidth=min(xwidth/(sol(5)*24.0),0.25)
      
      redj=0.09205 !radius of Earth divided by Jupiter
       
      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi   
     

      eccn=sqrt(sol(14)*sol(14)+sol(15)*sol(15)) !eccentricity

C       The case of the circular orbit is trivial
        tcentre=sol(7)/sol(5)-floor(sol(7)/sol(5))
        if(tcentre.lt.0.0d0) tcentre=tcentre+1.0d0!make sure between 0-1
        if(tcentre.gt.1.0d0) tcentre=tcentre-1.0d0
        ecentre=sol(7)/sol(5)-floor(sol(7)/sol(5))+0.5
        if(ecentre.lt.0.0d0) ecentre=ecentre+1.0d0!make sure between 0-1
        if(ecentre.gt.1.0d0) ecentre=ecentre-1.0d0

c      write(0,*) "tc,ec:",tcentre,ecentre

C     convert escale to real (used for eclipse scale)
      rescale=real(escale)
      
      call phasept(npt,time,phase,sol(5),0.0d0)
      do 5 i=1,npt
        phase(i)=phase(i)-tcentre
        if(phase(i).lt.-0.5d0) phase(i)=phase(i)+1.0d0
        if(phase(i).gt. 0.5d0) phase(i)=phase(i)-1.0d0
 5    continue
 
c      if(bins.gt.0) call phaseptsp(npt,time,phase,sol(5),0.0d0,sym)
      do 8 i=1,npt !copy data, because binning is destructive
        phase2(i)=phase(i)
        flux2(i)=flux(i)
        ferr2(i)=ferr(i)
 8    continue
      npt2=npt  
      write(0,*) "hello1"
      if(bins.gt.0) call binp(npt2,phase2,flux2,ferr2,bins,flag)!binning
      
      write(0,*) "hello2"
      
C     Find means of dataset
      itmean=0.0d0 !itmean is the mean exposure time
      bb(1)=time(1) !lower x-axis bound 
      bb(2)=time(1) !high x-axis bound
      do 11 i=1,npt
        bb(1)=min(time(i),bb(1)) !mix min and max times
        bb(2)=max(time(i),bb(2))
        itmean=itmean+itime(i)
 11   continue
      itmean=itmean/dble(npt) !mean exposure time
C     time of first observation - used to get the epoch below
      mintime=bb(1) 
 
C     Find bounds of data that will be plotted.
      bb(3)= 99.9d30!flux(1) !lower y-axis bound
      bb(4)=-99.9d30!flux(1) !high y-axis bound
      eoff2=-99.9d30!find max during secondary
      do 9 i=1,npt2
        if(abs(phase2(i)).lt.xwidth)then
            bb(3)=min(flux2(i),bb(3))
            bb(4)=max(flux2(i),bb(4))
        endif
        if((phase2(i).gt.0.5-xwidth).or.(phase2(i).lt.-0.5+xwidth))then
            eoff2=max(eoff2,flux2(i)-1.0)
        endif        
 9    continue
      if(eoff.le.0.0) eoff=2.0*abs(1.0-bb(4)) !auto offsets
c      eoff=2.0*abs(1.0-bb(4))  !remove this line
c      eoff=2.0*abs(1.0-bb(4))
      if(eoff2.le.0.0) eoff2=0.0
C     Choose between next two line for top of plot blank space
c      eoff2=eoff2+eoff2*0.15 !tag on 15% for nice plot boundaries
      eoff2=eoff
      if(nsensor.eq.2) then !when labels go away, more space for play
        bb(3)=bb(3)-(1.0-bb(3))*0.05
      else
        bb(3)=bb(3)-(1.0-bb(3))*0.25
      endif
      
      dnpt=dble(nplot-1) !convert integer to real for division 
      do 10 i=1,nplot  !generate 1 phase of model data
C     Psec contains the orbital Period of the planet in seconds
        ptime(i)=bb(1)+(bb(2)-bb(1))*dble(i-1)/dnpt !observation times
        if(bins.gt.0)then
            pitime(i)=sol(5)/dble(bins)
        else
            pitime(i)=itmean !mean integration times
        endif
        dtype(i)=0
 10   continue
            
      call pgpanl(npx,npy)
      bb(1)=-real(xwidth)
      bb(2)= real(xwidth)
      rbb(1)=real(bb(1))
      rbb(2)=real(bb(2))
      rbb(3)=real(bb(3)-0.35*(bb(4)-bb(3)))
      rbb(4)=real(eoff+eoff2)+1.0
      
      call pgwindow(24.0*rbb(1)*real(sol(5)),24.0*rbb(2)*real(sol(5)),
     .  rbb(3),rbb(4))
      
      call pgsch(1.5) !2.0
      if(nsensor.le.2) call pgslw(3)
      call pgbox('BCNTS1',0.0,0,'BNTSV1',0.0,0)
c      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgsch(2.0) !c1 made bigger labels.
      write(label,503) "Phase (Hours)"
 503  format(A13)
      call pgptxt((rbb(1)+rbb(2))/2.0,
     .  rbb(3)-0.14*(rbb(4)-rbb(3)),
     .  0.0,0.5,label)
      call pgptxt(rbb(1)-0.11*(rbb(2)-rbb(1)),(rbb(4)+rbb(3))/2.0,
     .  90.0,0.5,"Relative Flux")
      call pgslw(1)
      call pgsch(1.0)
      
      call phaseptsp(npt,time,phase,sol(5),0.0d0,sym)      
      if(bins.gt.0) call binp(npt,phase,flux,ferr,bins,flag)

      np=0      
      do 14 i=1,npt
        temp=phase(i)-tcentre
        if(temp.lt.-0.5) temp=temp+1.0
        if(temp.gt. 0.5) temp=temp-1.0
        if((temp.ge.bb(1)).and.
     .    (temp.le.bb(2)))then
            np=np+1
            px(np)=real(temp)!real(phase(i))-tcentre
            py(np)=real(flux(i))
            pz(np)=real(ferr(i))
            ps(np)=sym(i)
c            write(6,*) abs(px(np)),py(np),pz(np)
        endif
 14   continue
      call pgsch(2.0)
c      call pgerrb(6,np,px,py,pz,1.0)
      call pgpnts(np,px,py,ps,np)
      call pgsch(1.0)
      
      zpt=real(10.0**(sol(8)/-2.5))
      bscale=(rbb(4)-rbb(3))/(2.0*rescale)
      offset=((zpt+eoff)-(rbb(4)+rbb(3))/2.0)/rescale
      call pgwindow(rbb(1),rbb(2),zpt-bscale,zpt+bscale)
      call pgsch(2.0) !2.0
      if(nsensor.le.2) call pgslw(3)
      call pgptxt(rbb(2)+0.12*(rbb(2)-rbb(1)),zpt,
     .  90.0,0.5,"Folded Amplitude (ppm)")
      call pgslw(1)
      call pgsch(1.0)
c      call pgbox('',0.0,0,'CMTSV1',0.0,0)
      np=0      
      do 15 i=1,npt
        temp=phase(i)-ecentre
        if(temp.lt.-0.5) temp=temp+1.0
        if(temp.gt. 0.5) temp=temp-1.0
        if((temp.ge.bb(1)).and.
     .    temp.le.bb(2))then
            np=np+1
            px(np)=real(temp)!real(phase(i))-ecentre
            py(np)=real(flux(i))+offset
            pz(np)=real(ferr(i))
        endif
 15   continue
c      call pgerrb(6,np,px,py,pz,1.0)
      call pgpt(np,px,py,4)
c      write(6,*) 1.0-(zpt-bscale-offset),1.0-(zpt+bscale-offset)
      call pgwindow(rbb(1),rbb(2),1.0e6*(1.0-(zpt-bscale-offset)),
     .  1.0e6*(1.0-(zpt+bscale-offset)))
      call pgsch(1.5)!2.0
      if(nsensor.le.2) call pgslw(3)
      call pgbox('',0.0,0,'CMTSV1',0.0,0)
c1      call pgbox('',0.0,0,'CTSV1',0.0,0)
      call pgsch(1.0)
      call pgslw(1)
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))

      call transitmodel(nplot,ptime,pitime,dtype,tmodel,nfit,sol)
      call phasept(nplot,ptime,phase,sol(5),0.0d0)
      do 16 i=1,nplot
        terr(i)=1.0
 16   continue
 
      nplot2=nplot
c      call binp(nplot2,phase,tmodel,terr,bins,flag)

      do 6 i=1,nplot2
        phase(i)=phase(i)-tcentre
        if(phase(i).lt.-0.5) phase(i)=phase(i)+1.0
        if(phase(i).gt. 0.5) phase(i)=phase(i)-1.0
 6    continue
      call sort2(nplot2,phase,tmodel) !sort in phase
      
      np=0      
      do 12 i=1,nplot2
        temp=phase(i)
        if((temp.ge.bb(1)).and.
     .    (temp.le.bb(2)))then
            np=np+1
            px(np)=real(temp)!real(phase(i))-tcentre
            py(np)=real(10.0**(tmodel(i)/-2.5))
        endif
 12   continue

      call pgsch(2.0)
      call pgslw(3) !thicker lines
      call pgsci(2)
      call pgline(np,px,py)
      
      call pgwindow(rbb(1),rbb(2),zpt-bscale,zpt+bscale)
      
      do 7 i=1,nplot2
        phase(i)=phase(i)+tcentre-ecentre
        if(phase(i).lt.-0.5) phase(i)=phase(i)+1.0
        if(phase(i).gt. 0.5) phase(i)=phase(i)-1.0
 7    continue
      call sort2(nplot2,phase,tmodel) !sort in phase
      
      np=0      
      do 13 i=1,nplot2
        if((phase(i).ge.bb(1)).and.
     .    (phase(i).le.bb(2)))then
            np=np+1
            px(np)=real(phase(i))
            py(np)=real(10.0**(tmodel(i)/-2.5))+offset
        endif
 13   continue
      call pgsci(3)
      call pgline(np,px,py)      
      
      call pgsci(1.0)
      call pgslw(1)
      call pgsch(1.0)
 
C     Add labels.
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))      
      call pgsch(1.5)
      if(nsensor.le.2) call pgslw(3) !thicker lines
      write(label,500) "incl=",sol(6)
 500  format(A5,F5.2)
      write(6,*) "NSENSOR",nsensor
      if(nsensor.lt.2) call pgptxt(rbb(2)-0.03*(rbb(2)-rbb(1)),
     .  rbb(3)+0.055*(rbb(4)-rbb(3)),
     .  0.0,1.0,label)
      if(sol(4)/redj.lt.3.0d0)then
        write(label,501) "Rp=",sol(4)/redj,"Re"
      else
        write(label,501) "Rp=",sol(4),"Rj"
      endif
 501  format(A3,F5.2,1X,A2)
      if(nsensor.lt.2) call pgptxt(rbb(2)-0.03*(rbb(2)-rbb(1)),
     .  rbb(3)+0.11*(rbb(4)-rbb(3)),
     .  0.0,1.0,label)
      if(err(16).gt.0.0d0)then
        if(sol(16)/err(16).lt.10.0)then
            write(label,505) "Eclip-sig ",sol(16)/err(16)
        else
            write(label,512) "Eclip-sig ",sol(16)/err(16)
        endif
 505    format(A10,F4.1)
 512    format(A10,F5.1)
        if(nsensor.lt.2) call pgptxt(rbb(2)-0.03*(rbb(2)-rbb(1)),
     .      rbb(3)+0.275*(rbb(4)-rbb(3)),
     .      0.0,1.0,label)
      endif
      write(6,*) "doe:",doe
      if(doe.gt.0.0d0)then
        write(label,505) "Depth-sig ",doe
        if(nsensor.lt.2) call pgptxt(rbb(2)-0.03*(rbb(2)-rbb(1)),
     .      rbb(3)+0.220*(rbb(4)-rbb(3)),
     .      0.0,1.0,label)
      endif


      poff=0.5-sol(7)/tpi     

      epoch=sol(7)
      write(label,502) "Tc(E)=",epoch,"(BJD)+E(",sol(5),"days)"
 502  format(A6,F8.4,A8,F9.5,A5)
      if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .  rbb(3)+0.11*(rbb(4)-rbb(3)),
     .  0.0,0.0,label)
      write(label,504) "       ",err(7),"         ",err(5),"     "
 504  format(A7,F8.4,A9,F9.5,A5)
      if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .  rbb(3)+0.055*(rbb(4)-rbb(3)),
     .  0.0,0.0,label)
      
      if(sol(16)/err(16).ge.2.0d0)then !only plot if a reasonable detect
        if(Thpars(1,1).lt.10.0)then
            write(label,506) "Ag = ",Thpars(1,1),"(",Thpars(1,2),")"
 506        format(A5,F4.2,A1,F4.2,A1)
        else
            write(label,509) "Ag =",Thpars(1,1),"(",Thpars(1,2),")"
 509        format(A4,F6.2,A1,F6.2,A1)
        endif
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .      rbb(3)+0.275*(rbb(4)-rbb(3)),
     .      0.0,0.0,label)
     
        if(Thpars(2,1).lt.1.0d4)then
            write(label,507) "Teq=",int(Thpars(2,1)+0.5),"(",
     .          int(Thpars(2,2)+0.5),")"
 507        format(A4,I4,A1,I4,A1)
        else
            write(label,510) "Teq=",int(Thpars(2,1)+0.5),"(",
     .          int(Thpars(2,2)+0.5),")"
 510        format(A4,I5,A1,I5,A1)
        endif
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .      rbb(3)+0.220*(rbb(4)-rbb(3)),
     .      0.0,0.0,label)
        if(Thpars(3,1).lt.1.0d4)then
            write(label,508) "Teff=",int(Thpars(3,1)+0.5),"(",
     .          int(Thpars(3,2)+0.5),")"
 508        format(A5,I4,A1,I4,A1)
        else
            write(label,511) "Teff=",int(Thpars(3,1)+0.5),"(",
     .          int(Thpars(3,2)+0.5),")"
 511        format(A5,I5,A1,I5,A1)
        endif
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .      rbb(3)+0.165*(rbb(4)-rbb(3)),
     .      0.0,0.0,label)
      endif  
        
      call pgsch(1.0) !return character and line widths to normal
      call pgslw(1)   
        
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine panelA(npt,px,py,pz,time,flux,ferr,npx,npy,id,id2,kID,
     .  Kmag,Teff,logg,rad,Zerotime,nsensor,nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Used in transitplot to make TCERT plots
      implicit none
      integer npt,i,npx,npy,id,kID,Teff,id2,nsensor,nfit
      real px(npt),py(npt),bb(4),pz(npt),lx(2),ly(2)
      double precision time(npt),flux(npt),Kmag,logg,rad,Zerotime,
     .  ferr(npt),sol(nfit),Tmax,Tmin,T0,ep,p1
      character*80 label 
      
      call pgpanl(npx,npy)
      
      bb(1)=real(time(1))
      bb(2)=real(time(1))
      bb(3)=real(flux(1))
      bb(4)=real(flux(1))
      Tmax=time(1)
      Tmin=time(1)
      do 10 i=1,npt
        px(i)=real(time(i))
        py(i)=real(flux(i))
        pz(i)=real(ferr(i))
        bb(1)=min(bb(1),px(i))
        bb(2)=max(bb(2),px(i))
        bb(3)=min(bb(3),py(i))
        bb(4)=max(bb(4),py(i))
        Tmax=max(Tmax,time(i))
        Tmin=min(Tmin,time(i))
 10   continue
      bb(1)=floor(bb(1))
      bb(2)=floor(bb(2))+1.0
      bb(3)=bb(3)-0.1*(bb(4)-bb(3))
      bb(4)=bb(4)+0.1*(bb(4)-bb(3))
c      bb(1)=210.0
c      bb(2)=230.0
      
      call pgwindow(bb(1),bb(2),bb(3),bb(4))
C     range box, ticks and scale
      call pgsch(1.5) !2.0
      if(nsensor.le.2) call pgslw(3) !thicker lines
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
      write(label,502) "Time BJD-",int(Zerotime+2400000)
      write(6,*) "Zero: ",Zerotime
c 502  format(A9,F13.6)
c 502  format(A9,F16.8)
 502  format(A9,I7)
      call pgptxt(real((bb(1)+bb(2))/2.0),
     .  real(bb(3)-0.14*(bb(4)-bb(3))),
     .  0.0,0.5,label)
      call pgptxt(real(bb(1)-0.11*(bb(2)-bb(1))),real((bb(4)+bb(3))/2),
     .  90.0,0.5,"Relative Intensity")
      write(label,500) "Kepler:",id,".0",id2
c      write(label,500) "  MOST:",id,".0",id2     
 500  format(A7,I5,A2,I1)
      if(nsensor.lt.2) call pgptxt(real(bb(1)+(bb(2)-bb(1))/5.0),
     .  real(bb(4)+0.05*(bb(4)-bb(3))),
     .  0.0,0.5,label)
      if(nsensor.eq.0)then
        write(label,501) kID,Kmag,Teff,logg,rad
 501    format(I8,1X,F6.3,1X,I4,2(1X,F5.2))
c        write(label,504) Kmag,Teff,logg,rad
c 504    format(9X,F6.3,1X,I4,2(1X,F5.2))
      else
        write(label,503) floor(Kmag*10.0+0.5)/10.0,
     .      ((Teff+50)/100)*100,
     .      floor(logg*10.0+0.5)/10.0,
     .      floor(rad*10.0+0.5)/10.0
 503    format(8X,1X,F4.1,3X,I4,2(2X,F4.1))
      endif
      if(nsensor.lt.2) call pgptxt(real(bb(1)+2.0*(bb(2)-bb(1))/3.0),
     .  real(bb(4)+0.05*(bb(4)-bb(3))),
     .  0.0,0.5,label)
      
      call pgsci(2)
      ly(1)=1.000!bb(3)
      ly(2)=bb(4)
      T0=sol(7)-int((sol(7)-Tmin)/sol(5))*sol(5)
      write(6,*) "T0: ",T0
      do 11 i=1,int((Tmax-T0)/sol(5))+1
        lx(1)=real(T0)+real(i-1)*real(sol(5))
        lx(2)=lx(1)
        call pgsls(4)
        call pgline(2,lx,ly)
        call pgsls(1)
 11   continue
      call pgsci(1)
      
c      call pgsci(3)
c      p1=6.3111536272E+00
c      ep=1.0705818018E+02
c      T0=ep-int((ep-Tmin)/p1)*p1
c      write(6,*) "T0: ",T0
c      do 12 i=1,int((Tmax-T0)/p1)+1
c        lx(1)=real(T0)+real(i-1)*real(p1)
c        lx(2)=lx(1)
c        call pgsls(4)
c        call pgline(2,lx,ly)
c        call pgsls(1)
c 12   continue
c      call pgsci(1)
      
c      call pgsci(4)
c      p1=1.3478355224E+01
c      ep=6.9633837900E+01
c      T0=ep-int((ep-Tmin)/p1)*p1
c      write(6,*) "T0: ",T0
c      do 13 i=1,int((Tmax-T0)/p1)+1
c        lx(1)=real(T0)+real(i-1)*real(p1)
c        lx(2)=lx(1)
c        call pgsls(4)
c        call pgline(2,lx,ly)
c        call pgsls(1)
c 13   continue
c      call pgsci(1)
      
      call pgslw(1)
      call pgsch(0.8) !reduce marker size
      call pgpt(npt,px,py,17)
c      call pgerrb(6,npt,px,py,pz,1.0)
      call pgsci(1.0)
      
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine phaseptsp(npt,time,phase,period,toff,sym)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     npt - number of data points
C     time - times of observations
C     mag - the observations
C     phase - returned phase of observation
C     period - fixed period for data
      implicit none
      integer npt,sym(npt)
      double precision time(npt),phase(npt),period,toff

      integer i
      double precision temp

      do 10 i=1,npt
         temp=int(time(i)/period)
C        Get the phase
         phase(i)=time(i)/period-temp
         if(mod(temp+1,2.0).eq.0)then
            sym(i)=2
         else
            sym(i)=3
         endif
C        apply optional phase offset to make plot pretty
         phase(i)=phase(i)+toff
C        make sure phase is between 0 and 1
         if(phase(i).lt.0.0) phase(i)=phase(i)+1.0
         if(phase(i).gt.1.0) phase(i)=phase(i)-1.0
 10   continue
      return
      end



 
