CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine transitmodel(nfit,sol,npt,time,itime,tmodel,dtype)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,npt,i,j,nintg,dtype(npt)
      parameter(nintg=11)
      double precision sol(nfit),per,epoch,b,RpRs,tmodel(npt),
     .  time(npt),phi,adrs,bt(nintg),bs2,Pi,tpi,c1,c2,c3,c4,
     .  itime(npt),t,tflux(nintg),dnintg,dnintgm1,pid2,zpt,eccn,w,Eanom,
     .  Tanom,trueanomaly,phi0,vt(nintg),K,voff,drs,distance,dil,ab,
     .  ratio,ted,fDB,Cs,tide(nintg),ell,Manom,ltt,alb(nintg),ag,
     .  albedomod,phase
      
      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi 
      pid2=Pi/2.0d0
      Cs=2.99792458e8 !Speed of light
      fDB=1.896 !Doppler Boosting factor
      
      ltt=39.1709d0/(86400.0d0) !light travel time
      
      ell=sol(16)/1.0d6 !Ellipsoidal variations
      ag=sol(17)/1.0d6 !Phase changes
      
      per=sol(2)     !Period (days)
      b=sqrt(sol(3))       !impact parameter
      bs2=sol(3)
      RpRs=sol(4)    !Rp/R*
      zpt=sol(6)
      eccn=sqrt(sol(7)*sol(7)+sol(8)*sol(8)) !eccentricity
      if(eccn.ge.1.0) eccn=0.99
      if(eccn.eq.0.0d0)then
        w=0.0d0
      else
        w=atan(sol(8)/sol(7))
        if((sol(7).gt.0.0d0).and.(sol(8).lt.0.0d0))then
            w=tPi+w
        elseif((sol(7).lt.0.0d0).and.(sol(8).ge.0.0d0))then 
            w=Pi+w
        elseif((sol(7).lt.0.0d0).and.(sol(8).lt.0.0d0))then
            w=Pi+w
        endif
      endif
c      write(0,*) sol(7),sol(8),w
c      write(0,*) "w:",acos(sol(7)/eccn),asin(sol(8)/eccn)
c      read(5,*)

C     Depth of secondary
      ted=sol(15)/1.0d6

C     a/R*
c      adrs=sol(5)*per/tpi*sqrt(1-sol(3))*(1+sol(8))/sqrt(1-eccn*eccn)
      adrs=sol(5)*per/tpi*sqrt((1.0d0+sol(4))**2.0d0-sol(3))*
     .  (1+sol(8))/sqrt(1-eccn*eccn)

      K=sol(9)
      voff=sol(10)
      c1=sol(11)      !non-linear limb-darkening
      c2=sol(12)
      c3=sol(13)
      c4=sol(14)
      
      dil=sol(18) !dilution parameter (model scaling)
      
      dnintg=dble(nintg) !convert integer to double
      dnintgm1=2.0*dnintg-2.0
      
C     Find phase at centre of transit
c      phi0=tpi*(epoch/per-floor(epoch/per))
c      call kepler(phi0,Eanom,eccn)
c      Tanom=trueanomaly(eccn,Eanom)
c      phi0=Tanom+w
      Eanom=w
      epoch=sol(1)   !center of transit time (days)
      Manom=w !mean anomaly
      call kepler(Manom,Eanom,eccn) !eccentric anomaly
      phi0=trueanomaly(eccn,Eanom)
      
      do 10 i=1,npt
        do 11 j=1,nintg
c            itime(i)=0.0
            tflux(j)=0.0 !initialize model
C           sample over integration time
            t=time(i)+itime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1-epoch
C           get orbital position (mean anomaly)
c            phi=(t-epoch)/per-floor((t-epoch)/per)
            phi=t/per-floor(t/per)
            phi=phi*tPi
c            call kepler(phi,Eanom,eccn) !get eccentric anomaly
c            Tanom=trueanomaly(eccn,Eanom) !get true anonaly
            Manom=phi+w
            if(Manom.gt.tPi) Manom=Manom-tPi
            if(Manom.lt.0.0d0) Manom=Manom+tPi
            call kepler(Manom,Eanom,eccn)
            Tanom=trueanomaly(eccn,Eanom)
            if(phi.gt.Pi) phi=phi-tPi            
            drs=distance(adrs,eccn,Tanom)
            bt(j)=sqrt(bs2+(drs*sin(Tanom-phi0))**2)
C           Correct for light-travel time!
c            if((abs(bt(j))-RpRs.le.1.0d0).and.(abs(phi).gt.Pid2))then
c                t=time(i)-ltt+itime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1
c     .              -epoch
c                phi=(t-ltt)/per-floor((t-ltt)/per)
c                phi=phi*tPi
c                Manom=phi+w
c                if(Manom.gt.tPi) Manom=Manom-tPi
c                if(Manom.lt.0.0d0) Manom=Manom+tPi
c                call kepler(Manom,Eanom,eccn)
c                Tanom=trueanomaly(eccn,Eanom)
c                if(phi.gt.Pi) phi=phi-tPi            
c                drs=distance(adrs,eccn,Tanom)
c                bt(j)=sqrt(bs2+(drs*sin(Tanom-phi0))**2)
c            endif
            vt(j)=K*(cos(Pid2+Tanom-phi0)+eccn*cos(w))
            tide(j)=ell*cos(2.0d0*(Pid2+Tanom-phi0))
            alb(j)=albedomod(Pi,ag,Tanom-phi0)
            
            if(j.eq.nintg/2+1)then
                phase=Tanom-phi0!phi(nintg/2+1)
                if(phase.gt.Pi) phase=phase-tPi
                if(phase.lt.-Pi) phase=phase+tPi
            endif
c            write(0,*) tide(j),ell,cos(2.0d0*(Tanom-phi0))
c            read(5,*)
c            vt(j)=K*(cos(Pid2+w+Tanom-phi0)+eccn*cos(Pid2+Tanom))
c            write(6,*) t,time(i),itime(i)
c            read(5,*)
 11     continue
        if(dtype(i).eq.0)then
C           compute mid-exposure orbital position (eclipse of transit?)
c            phi=(time(i)-epoch)/per-floor((time(i)-epoch)/per)
c            phi=phi*tPi
c            if(phi.gt.Pi) phi=phi-tPi
c            if(abs(phi).lt.Pid2)then
            if(abs(phase).lt.Pid2)then
C       If we have a transit
                call occultsmall(RpRs,c1,c2,c3,c4,nintg,bt,tflux)
                tmodel(i)=0.0d0
                do 12 j=1,nintg
C                   model=transit+doppler+ellipsodial 
                    tmodel(i)=tmodel(i)+tflux(j)-fDB*vt(j)/Cs+tide(j)+
     .                  alb(j)
 12             continue
                tmodel(i)=tmodel(i)/dnintg
            else
C       We have an eclipse
                tmodel(i)=0.0d0
                do 14 j=1,nintg
                    ratio=1.0d0
                    ab=dabs(bt(j))
                    if((ab.ge.1.0d0).and.(ab-RpRs.le.1.0d0))then
                        ratio=(1.0d0+RpRs-ab)/(2.0d0*RpRs)
                    elseif((ab.lt.1.0d0).and.(ab+RpRs.ge.1.0d0))then
                        ratio=(RpRs+1.0d0-ab)/(2.0d0*RpRs)
                    elseif(ab-RpRs.gt.1.0d0)then
                        ratio=0.0d0
                    endif
c                    write(0,*) ab,RpRs,ratio
c                    read(5,*) 
                    tmodel(i)=tmodel(i)+(1.0d0-ted*ratio)
     .                  -fDB*vt(j)/Cs+tide(j)+alb(j)
 14             continue
                tmodel(i)=tmodel(i)/dnintg
            endif
            tmodel(i)=tmodel(i)+(1.0-tmodel(i))*dil+zpt !add dilution
         else
            tmodel(i)=0.0d0
            do 13 j=1,nintg
                tmodel(i)=tmodel(i)+vt(j)
 13         continue
            tmodel(i)=tmodel(i)/dnintg+voff
c            write(0,*) "rv:",tmodel(i)
c            read(5,*)
         endif
 10   continue
      
c      do 15 i=1,npt
c        write(6,*) time(i),tmodel(i)
c 15   continue
      
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