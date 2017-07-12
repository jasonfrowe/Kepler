CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine transitmodel(nfit,sol,npt,time,itime,tmodel)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,npt,i,j,nintg
      parameter(nintg=11)
      double precision sol(nfit),per,epoch,b,RpRs,tmodel(npt),
     .  time(npt),phi,adrs,bt(nintg),bs2,Pi,tpi,c1,c2,c3,c4,
     .  itime(npt),t,tflux(nintg),dnintg,dnintgm1,pid2,zpt,eccn,w
      
      per=sol(2)     !Period (days)
      b=sol(3)       !impact parameter
      RpRs=sol(4)    !Rp/R*
      epoch=sol(1)   !center of transit time (days)
      adrs=sol(5)    !a/R*
      zpt=sol(6)
      eccn=sqrt(sol(7)*sol(7)+sol(8)*sol(8)) !eccentricity
      if(eccn.eq.0.0d0)then
        w=0.0d0
      else
        w=acos(sol(7)/eccn)
      endif
      c1=sol(9)      !non-linear limb-darkening
      c2=sol(10)
      c3=sol(11)
      c4=sol(12)
      
      dnintg=dble(nintg) !convert integer to double
      dnintgm1=2.0*dnintg-2.0
      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi 
      pid2=Pi/2.0d0
      bs2=b*b
      
      do 10 i=1,npt
        do 11 j=1,nintg
            tflux(j)=0.0 !initialize model
C           sample over integration time
            t=time(i)+itime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1
C           get orbital position
            phi=(t-epoch)/per-floor((t-epoch)/per)
            phi=phi*tPi
            if(phi.gt.Pi) phi=phi-tPi
C           determine spatial impact parameter
            bt(j)=sqrt(bs2+(adrs*sin(phi))**2)
c            write(6,*) t,time(i),itime(i)
c            read(5,*)
 11     continue
C       compute mid-exposure orbital position (eclipse of transit?)
        phi=(time(i)-epoch)/per-floor((time(i)-epoch)/per)
        phi=phi*tPi
        if(phi.gt.Pi) phi=phi-tPi
        if(abs(phi).lt.Pid2)then
C       If we have a transit
            call occultsmall(RpRs,c1,c2,c3,c4,nintg,bt,tflux)
            tmodel(i)=0.0d0
            do 12 j=1,nintg
                tmodel(i)=tmodel(i)+tflux(j)
 12         continue
            tmodel(i)=tmodel(i)/dnintg+zpt
        else
C       We have an eclipse
            tmodel(i)=1.0d0+zpt
        endif
 10   continue
      
c      do 13 i=1,npt
c        write(6,*) time(i),tmodel(i)
c 13   continue
      
      
      return
      end