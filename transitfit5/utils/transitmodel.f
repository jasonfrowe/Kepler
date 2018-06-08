CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,
     .  itime,ntt,tobs,omc,tmodel,dtype)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     My Transit-model.  Handles multi-planets, TTVs, phase-curves and
C     radial-velocities
C     (c) jasonfrowe@gmail.com
      implicit none
      integer nfit,npt,i,j,nintg,dtype(npt),ii,nplanet,nplanetmax,nmax,
     .   caltran
      parameter(nintg=41)
      double precision sol(nfit),per,epoch,b,RpRs,tmodel(npt),
     .  time(npt),phi,adrs,bt(nintg),bs2,Pi,tpi,c1,c2,c3,c4,
     .  itime(npt),t,tflux(nintg),dnintg,dnintgm1,pid2,zpt,eccn,w,Eanom,
     .  Tanom,trueanomaly,phi0,vt(nintg),K,voff,drs,distance,dil,G,
     .  esw,ecw,Manom,Ted,ell,Ag,Cs,fDB,tide(nintg),alb(nintg),
     .  albedomod,phase,ratio,ab,tm,mu(nintg),y2,x2,incl,mulimbf(nintg),
     .  occ(nintg),bp(nintg),jm1,tdnintg
      integer ntt(nplanetmax)
      double precision tobs(nplanetmax,nmax),omc(nplanetmax,nmax),ttcor
      
      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi 
      pid2=Pi/2.0d0
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      Cs=2.99792458e8 !Speed of light
c      fDB=2.21 !Doppler Boosting factor
      fDB=1.0
      
      c1=sol(2)      !non-linear limb-darkening
      c2=sol(3)
      c3=sol(4)
      c4=sol(5)
      dil=sol(6)     !dilution parameter (model scaling)
      voff=sol(7)    !velocity zero point
      zpt=sol(8)     !flux zero point.
      
      do 17 i=1,npt
        tmodel(i)=0.0d0
 17   continue
      
      do 16 ii=1,nplanet
      
        per=sol(10*(ii-1)+8+2)     !Period (days)
c        sol(10*(ii-1)+8+3)=abs(sol(10*(ii-1)+8+3)) !make b positive
c        bs2=abs(sol(10*(ii-1)+8+3))
        bs2=sol(10*(ii-1)+8+3)*sol(10*(ii-1)+8+3)
        b=sqrt(bs2)       !impact parameter
c        write(0,*) bs2,b
        RpRs=abs(sol(10*(ii-1)+8+4))    !Rp/R*

        ecw=sol(10*(ii-1)+8+6)
        esw=sol(10*(ii-1)+8+5)
c        eccn=sqrt(ecw*ecw+esw*esw) !eccentricity
        eccn=(ecw*ecw+esw*esw)

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
      
c        write(0,*) sol(7),sol(8),w
c        write(0,*) "w:",acos(sol(7)/eccn),asin(sol(8)/eccn)
c        read(5,*)

C       a/R*
c        adrs=sol(5)*per/tpi*sqrt(1-sol(3))*(1+sol(8))/sqrt(1-eccn*eccn)
        adrs=1000.0*sol(1)*G*(Per*86400.0d0)**2/(3.0d0*Pi)
        adrs=adrs**(1.0d0/3.0d0)
c        write(0,*) "a/R*:",adrs

C       Find inclination !hmm.. this is probably wrong
        incl=acos(b/adrs)


c        K=abs(sol(10*(ii-1)+8+7))
        K=sol(10*(ii-1)+8+7)
      
        ted=sol(10*(ii-1)+8+8)/1.0d6 !Occultation Depth
        ell=sol(10*(ii-1)+8+9)/1.0d6 !Ellipsoidal variations
        ag=sol(10*(ii-1)+8+10)/1.0d6 !Phase changes
      
        dnintg=dble(nintg) !convert integer to double
        tdnintg=2.0d0*dnintg
        dnintgm1=2.0*dnintg-2.0
      
C     Find phase at centre of transit
        epoch=sol(10*(ii-1)+8+1)   !center of transit time (days)
c        phi1=(epoch/per-int(epoch/per))*twopi
        Eanom=tan(w/2.0d0)/sqrt((1.0d0+eccn)/(1.0d0-eccn)) !mean anomaly
        Eanom=2.0d0*atan(Eanom)
        phi0=Eanom-eccn*sin(Eanom)
      

!Add parallel commands here
!$OMP PARALLEL DO PRIVATE(j,jm1,ttcor,tflux,t,phi,Manom,Tanom,drs,incl,
!$OMP& x2,y2,bt,vt,tide,alb,caltran,mu,tm,bp,ratio,occ) 
!$OMP& FIRSTPRIVATE (Eanom,c1,c2,c3,c4)
        do i=1,npt
            call lininterp(tobs,omc,nplanetmax,nmax,ii,ntt,time(i),
     .          ttcor)
c            write(0,*) ii,time(i),ttcor
c            read(5,*)
            do 11 j=1,nintg
                jm1=dble(j-1)
                tflux(j)=0.0 !initialize model
C               sample over integration time

c                t=t-ttcor

!              old time-convolution
c                t=time(i)+itime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1-
c     .              epoch-ttcor

!              new time-convolution (basically gives same results)
                t=time(i)-itime(i)*(0.5d0-1.0d0/tdnintg-jm1/dnintg)-
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
C              Added this (2014/04/23)
                incl=acos(b/drs)
                x2=drs*Sin(Tanom-w)
                y2=drs*Cos(Tanom-w)*cos(incl)

c                x2=drs*Cos(Tanom+w)
c                y2=drs*Sin(Tanom+w)*cos(incl)


c                bt(j)=sqrt(bs2+(drs*sin(Tanom-w))**2)
                bt(j)=sqrt(x2*x2+y2*y2)
C               Correct for light-travel time!
c            if((abs(bt(j))-RpRs.le.1.0d0).and.(abs(phi).gt.Pid2))then
c              t=time(i)-ltt+itime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1
c     .              -epoch
c              phi=(t-ltt)/per-floor((t-ltt)/per)
c              phi=phi*tPi
c              Manom=phi+w
c              if(Manom.gt.tPi) Manom=Manom-tPi
c              if(Manom.lt.0.0d0) Manom=Manom+tPi
c              call kepler(Manom,Eanom,eccn)
c              Tanom=trueanomaly(eccn,Eanom)
c                if(phi.gt.Pi) phi=phi-tPi            
c                drs=distance(adrs,eccn,Tanom)
c                bt(j)=sqrt(bs2+(drs*sin(Tanom-phi0))**2)
c            endif
c                vt(j)=K*(cos(Pid2+Tanom-phi0)+eccn*cos(w))

                vt(j)=K*(cos(Tanom-w+pid2)+eccn*cos(-w+pid2))
c                write(6,*) "DB:",fDB*vt(j)/Cs,vt(j)
c                vt(j)=-K*(cos(Tanom+w)+ecw)

                tide(j)=ell*(drs/adrs)**(1.0d0/3.0d0)*
     .              cos(2.0d0*(Pid2+Tanom-w))
c                tide(j)=ell*(drs/adrs)**(1.0d0/3.0d0)*
c     .              cos(2.0d0*(Pid2+phi))
     
                alb(j)=albedomod(Pi,ag,Tanom-w)*adrs/drs
c                alb(j)=albedomod(Pi,ag,phi)*adrs/drs
            
c                if(j.eq.nintg/2+1)then
c                    phase=Tanom-w!phi(nintg/2+1)
c                    if(phase.gt.Pi) phase=phase-tPi
c                    if(phase.lt.-Pi) phase=phase+tPi
c                    write(6,*) time(i),x2/adrs,y2/adrs/cos(incl)
c                    write(6,*) time(i),vt(j),w
c                     write(6,*) time(i),vt(j),w
c                endif

c               write(6,*) t,x2,y2

 11         continue
            if(dtype(i).eq.0)then
c                if(abs(phase).lt.Pid2)then
                if(y2.ge.0.0d0)then
C       If we have a transit
                    caltran=0 !if zero, there is no transit
                    do 18 j=1,nintg
                        if(bt(j).le.1.0d0+RpRs)then
                           caltran=1
                        endif
 18                 continue
                    if(caltran.eq.1) then
                       !quadratic co-efficients
                       if((c3.eq.0.0).and.(c4.eq.0.0))then
                         call occultquad(bt,c1,c2,RpRs,tflux,mu,nintg)
c                         write(0,*) (tflux(j),j=1,nintg)
c                         write(0,*) (bt(j),j=1,nintg)
c                         read(5,*)
                      !Kipping co-efficients
                       elseif((c1.eq.0.0).and.(c2.eq.0.0))then
                          c1=2.0d0*sqrt(c3)*c4 !convert to regular LD
                          c2=sqrt(c3)*(1.0d0-2.0d0*c4)
                          call occultquad(bt,c1,c2,RpRs,tflux,mu,nintg)
                          c1=0.0d0  !zero out entries.
                          c2=0.0d0
                       else
                      !non-linear law.
                       call occultsmall(RpRs,c1,c2,c3,c4,nintg,bt,tflux)
c                         call occultnl(RpRs,c1,c2,c3,c4,bt,tflux,
c     .                     mulimbf,nintg)
c                         write(0,*) RpRs,abs(sol(10*(ii-1)+8+4))
c                         write(6,550) RpRs,(bt(j),tflux(j),j=1,nintg)
 550                     format(30(F8.5,1X))
                       endif
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
        enddo
!$OMP END PARALLEL DO
 
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
