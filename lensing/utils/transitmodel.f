CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,
     .  itime,ntt,tobs,omc,tmodel,dtype)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,npt,i,j,nintg,dtype(npt),ii,nplanet,nplanetmax,nmax,
     .   caltran,nrad,ntheta
      parameter(nintg=11)
      double precision sol(nfit),per,epoch,b,RpRs,tmodel(npt),
     .  time(npt),phi,adrs,bt(nintg),bs2,Pi,tpi,c1,c2,c3,c4,
     .  itime(npt),t,tflux(nintg),dnintg,dnintgm1,pid2,zpt,eccn,w,Eanom,
     .  Tanom,trueanomaly,phi0,vt(nintg),K,voff,drs,distance,dil,G,
     .  esw,ecw,Manom,Ted,Ag,Cs,fDB,tide(nintg),alb(nintg),
     .  albedomod,phase,ratio,ab,tm,mu(nintg),incl,mulimbf(nintg),
     .  occ(nintg),bp(nintg),fEL,Rsun,Rearth,Msun,Mearth,Ms,Mp,Psec,
     .  asemi,aConst,AU,Rp,Rs,eps,F0,ZFlux,nl(4),x2(nintg),y2(nintg),
     .  z2(nintg),Rstar,Rl,Ml,thres,Pflux,rsdau
      integer ntt(nplanetmax)
      double precision tobs(nplanetmax,nmax),omc(nplanetmax,nmax),ttcor
      
      nrad=5000   !number of rings to integrate star/planet surface
      ntheta=5000 !maximum number of angular divisions
      thres=1.0e-7 !tells us when amplification is negliable

      Rsun=696265.0d0*1000.0d0 !m
      Rearth=6371.0*1000.0d0 !m
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi 
      pid2=Pi/2.0d0
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      Cs=2.99792458e8 !Speed of light
      AU=1.49598e11

      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

      fDB=sol(10) !Doppler Boosting factor
      fEL=sol(11) !Ellipsoidal amplitude factor
      
      c1=sol(3)      !non-linear limb-darkening
      c2=sol(4)
      c3=sol(5)
      c4=sol(6)
      nl(1)=sol(3)
      nl(2)=sol(4)
      nl(3)=sol(5)
      nl(4)=sol(6)
      dil=sol(7)     !dilution parameter (model scaling)
      voff=sol(8)    !velocity zero point
      zpt=sol(9)     !flux zero point.

      F0=Zflux(nrad,ntheta,nl,Pi) !need to call it here for nl change.
c      write(0,*) "F0: ",F0
      
      do 17 i=1,npt
        tmodel(i)=0.0d0
 17   continue
      
      do 16 ii=1,nplanet
      
        per=sol(9*(ii-1)+11+2)     !Period (days)
        Psec=Per*8.64d4 !sec ; period of planet
c        sol(10*(ii-1)+8+3)=abs(sol(10*(ii-1)+8+3)) !make b positive
c        bs2=abs(sol(10*(ii-1)+8+3))
        bs2=sol(9*(ii-1)+11+3)*sol(9*(ii-1)+11+3)
        b=sqrt(bs2)       !impact parameter
c        write(0,*) bs2,b
        Rp=Rearth*abs(sol(9*(ii-1)+11+5))
        Rs=Rsun*abs(sol(2))
        rsdau=Rs/AU
        Rstar=abs(sol(2))
        Rl=Rp/Rsun !radius of lens (Rsun)
        RpRs=Rp/Rs!Rp/R*
c        write(0,*) "RpRs: ",RpRs,sol(2),sol(9*(ii-1)+11+5)
c        read(5,*)

        ecw=sol(9*(ii-1)+11+7)
        esw=sol(9*(ii-1)+11+6)
        eccn=sqrt(ecw*ecw+esw*esw) !eccentricity
        if(eccn.ge.1.0) eccn=0.99
        if(eccn.eq.0.0d0)then
            w=0.0d0
        else
            if(ecw.eq.0.0d0)then
                w=Pi/2.0d0
            else
                w=atan(esw/ecw)
            endif
            if((ecw.gt.0.0d0).and.(esw.lt.0.0d0))then
                w=tPi+w
            elseif((ecw.lt.0.0d0).and.(esw.ge.0.0d0))then 
                w=Pi+w
            elseif((ecw.le.0.0d0).and.(esw.lt.0.0d0))then
                w=Pi+w
            endif
        endif     

C       a/R*
        Ms=abs(sol(1))*Msun !kg ; mass of star
        Mp=abs(sol(9*(ii-1)+11+4))*Mearth !kg ; mass of planet
        ml=Mp/Msun !mass of lens (Msun)
        asemi=(Ms+Mp)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst
        adrs=asemi/Rs
c        adrs=1000.0*0.2*G*(Per*86400.0d0)**2/(3.0d0*Pi)
c        adrs=adrs**(1.0d0/3.0d0)
c        write(0,*) "a/R*:",adrs,asemi/AU
c        read(5,*)

C       Find inclination
        incl=acos(b/adrs)

        K=2.0*pi*G*Mp**3*(sin(incl))**3/
     .    (Psec*(1.0d0-eccn*eccn)**(3.0d0/2.0d0)*(Ms+Mp)*(Ms+Mp))
        K=K**(1.0d0/3.0d0)
c        K=sol(10*(ii-1)+8+7)
c        write(0,*) "K:",K
c        read(5,*)
      
        ag=sol(9*(ii-1)+11+8)/1.0d6 !Phase changes
        ted=sol(9*(ii-1)+11+9)/1.0d6 !Occultation Depth
        eps=Mp/Ms*(Rs/asemi)**3.0d0
c        ell=sol(9*(ii-1)+8+9)/1.0d6 !Ellipsoidal variations

      
        dnintg=dble(nintg) !convert integer to double
        dnintgm1=2.0*dnintg-2.0
      
C     Find phase at centre of transit
        epoch=sol(9*(ii-1)+11+1)   !center of transit time (days)
c        phi1=(epoch/per-int(epoch/per))*twopi
        Eanom=tan(w/2.0d0)/sqrt((1.0d0+eccn)/(1.0d0-eccn)) !mean anomaly
        Eanom=2.0d0*atan(Eanom)
        phi0=Eanom-eccn*sin(Eanom)
      
        do 10 i=1,npt
            call lininterp(tobs,omc,nplanetmax,nmax,ii,ntt,time(i),
     .          ttcor)

            do 11 j=1,nintg
                tflux(j)=0.0 !initialize model
C               sample over integration time
                t=time(i)+itime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1-
     .              epoch-ttcor

     
c                write(0,*) itime(i)
C               get orbital position (mean anomaly)
                phi=t/per-floor(t/per)
                phi=phi*tPi+phi0
                Manom=phi
                if(Manom.gt.tPi) Manom=Manom-tPi
                if(Manom.lt.0.0d0) Manom=Manom+tPi
                call kepler(Manom,Eanom,eccn)
                Tanom=trueanomaly(eccn,Eanom)
                if(phi.gt.Pi) phi=phi-tPi            
                drs=distance(adrs,eccn,Tanom)
                x2(j)=drs*Sin(Tanom-w)
                y2(j)=drs*Cos(Tanom-w)*cos(incl)
                z2(j)=drs*Cos(Tanom-w)

                bt(j)=sqrt(x2(j)*x2(j)+y2(j)*y2(j))

                vt(j)=K*(cos(Tanom-w+pid2)+eccn*cos(-w+pid2))

                tide(j)=fEL*eps*(drs/adrs)**(1.0d0/3.0d0)*
     .              cos(2.0d0*(Pid2+Tanom-w))

                alb(j)=albedomod(Pi,ag,Tanom-w)*adrs/drs
            

 11         continue
            if(dtype(i).eq.0)then
c                if(abs(phase).lt.Pid2)then
                if(y2(1).ge.0.0d0)then
C       If we have a transit
                    caltran=0 !if zero, there is no transit
                    do 18 j=1,nintg
                        if(bt(j).lt.1.0d0+RpRs*2.0d0)then
                           caltran=1
                        endif
 18                 continue
                    if(caltran.eq.1) then
                       if((c3.eq.0.0).and.(c4.eq.0.0))then
                         call occultquad(bt,c1,c2,RpRs,tflux,mu,nintg)
                       else
                       call occultsmall(RpRs,c1,c2,c3,c4,nintg,bt,tflux)
 550                     format(30(F8.5,1X))
                       endif
c                        do 21 j=1,nintg
c                           tflux(j)=PFlux(nrad,ntheta,rsdau*x2(j),
c     .                        rsdau*y2(j),rsdau*z2(j),
c     .                        Rstar,Rl,Ml,nl,Pi,G,Cs,Msun,Rsun,AU,F0,
c     .                        thres)/F0
c                           write(0,*) time(i),x2(j),tflux(j)
c 21                     continue
                    else
                        do 19 j=1,nintg
                           tflux(j)=1.0d0
 19                     continue
                    endif
                    tm=0.0d0
                    do 12 j=1,nintg
                        if(RpRs.le.0.0)tflux(j)=1.0d0
C                   model=transit+doppler+ellipsodial 
                        tm=tm+tflux(j)-fDB*vt(j)/Cs+tide(j)+alb(j)
 12                 continue
                    tm=tm/dnintg
                else
C       We have an eclipse
                    tm=0.0d0
                    do 20 j=1,nintg
                      bp(j)=bt(j)/RpRs
 20                 continue
                    call occultuniform(bp,1.0/RpRs,occ,nintg)
                    do 14 j=1,nintg
                        ratio=1.0d0-occ(j)

C                      Old estimate, replaced by analytic function
c                        ratio=1.0d0
c                        ab=dabs(bt(j))
c                        if((ab.ge.1.0d0).and.(ab-RpRs.le.1.0d0))then
c                            ratio=(1.0d0+RpRs-ab)/(2.0d0*RpRs)
c                        elseif((ab.lt.1.0d0).and.(ab+RpRs.ge.1.0d0))then
c                            ratio=(RpRs+1.0d0-ab)/(2.0d0*RpRs)
c                        elseif(ab-RpRs.gt.1.0d0)then
c                            ratio=0.0d0
c                        endif
c                        write(0,*) bt(j),ratio,rationew
c                        read(5,*)
                        if(RpRs.le.0.0d0) ratio=0.0d0
                        tm=tm+(1.0d0-ted*ratio)
     .                      -fDB*vt(j)/Cs+tide(j)+alb(j)
 14                 continue
                    tm=tm/dnintg
                endif
                tm=tm+(1.0d0-tm)*dil-1.0d0!add dilution
            else
                tm=0.0d0
                do 13 j=1,nintg
                    tm=tm+vt(j)
 13             continue
                tm=tm/dnintg
c            write(0,*) "rv:",tmodel(i)
c            read(5,*)
            endif
            tmodel(i)=tmodel(i)+tm
 10     continue
 
C     Need to add zero points (voff, zpt)
      
c        do 9 i=1,npt
c            write(6,*) time(i),tmodel(i)
c 9      continue
  
 16   continue  
 
      do 15 i=1,npt
        if(dtype(i).eq.0)then
            tmodel(i)=tmodel(i)+zpt+1.0d0
        else
            tmodel(i)=tmodel(i)+voff
        endif
 15   continue
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function albedomod(Pi,ag,phi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision Pi,phi,alpha,phase,ag

      phi=phi+Pi
      if(phi.gt.2.0*Pi) phi=phi-2.0*Pi


      alpha=abs(phi)      
c      alpha=2.0*Pi*t/Per+phi
      alpha=alpha-2.0*Pi*int(alpha/(2.0*Pi))
      if(alpha.gt.Pi) alpha=abs(alpha-2.0*pi)
c      write(6,*) t,alpha
c      phase=(1.0d0+cos(alpha))/2.0d0
      phase=(sin(alpha)+(Pi-alpha)*cos(alpha))/Pi  !Lambertian Sphere
      
      albedomod=ag*phase
      
      return
      end
