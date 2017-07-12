CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine panel(n)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n
      
      if (n.eq.1) then
        call pgsvp(0.1,0.95,0.83,0.98)
      elseif(n.eq.2) then
        call pgsvp(0.1,0.95,0.60,0.75)
      elseif(n.eq.3) then
        call pgsvp(0.1,0.50,0.10,0.50)
      elseif(n.eq.4) then
        call pgsvp(0.60,0.95,0.10,0.50)
      endif
      
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine opengraphics()
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calls PGPLOT routines to open graphic window for plotting
      implicit none
      
      call pgopen('?')
c      call pgopen('/null')
c      call pgpap(6.0,1.0)
c      call pgask(.false.)
      call pgsch(0.8)
c      call pgsubp(1,4)
c      call pgvport(0.1,0.95,0.2,0.8) !gives enough room for labels
      
      call PGPAP(8.0,0.8)

c      call pgpage()
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine closegraphics()
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     close up graphics device
      implicit none
      
      call pgclos()
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotdata(npt,time,mag,merr,x1,x2,pan,MOSTtime,nmax,px,
     .  py,pz,bb)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine will plot the data (ie. time vs mag)
C     integer npt: Total number of data points
C     real*8 time(npt): The x-axis data points
C     real*8 mag(npt): The y-axis data points
C     real*8 merr(npt): Error associated with mag()
C     real*8 x1,x2: range of data to plot
C     integer panx,pany: PGPLOT panel to plot data in
C     integer nmax: maximum number of observations expected
C     real px(nmax),py(nmax) real arrays for plotting
C     real*8 bb(4) returns plotting bounds used
      implicit none
      integer npt,i,sym,nmax,pan
      real px(nmax),py(nmax),pz(nmax)
      double precision time(npt),mag(npt),merr(npt),x1,x2,mmin,mmax,
     .  tmin,tmax,bb(4),tmp(2)
      double precision MOSTtime
      character*80 xlabel
      
C     If we have few observations, use big dots, otherwise, small dots
      sym=17
      if(npt.gt.1000) sym=-1

C     select panel for plot
c      call pgpanl(panx,pany)
      call panel(pan)
      
C     find the boundaries of the data
      mmin=mag(1)
      mmax=mag(1)
      tmin=time(1)
      tmax=time(1)
      do 10 i=2,npt
         if(mag(i).gt.mmax) mmax=mag(i)
         if(mag(i).lt.mmin) mmin=mag(i)
         if(time(i).gt.tmax) tmax=time(i)
         if(time(i).lt.tmin) tmin=time(i)
 10   continue
C     Plots will be in mmag
      mmin=mmin*1000.0
      mmax=mmax*1000.0

      if(x1.ne.-1.0) then
         tmin=x1
      endif
      if(x2.ne.-1.0) then
         tmax=x2
      endif

      write(0,*) "---------------------"
      write(0,*) "tmin, tmax",tmin,tmax
      write(0,*) "mmin, mmax",mmin,mmax
      write(0,*) "---------------------"
C     Set plotting range
      tmp(1)=mmin-0.1*(mmax-mmin)
      tmp(2)=mmax+0.1*(mmax-mmin)
      mmin=tmp(1)
      mmax=tmp(2)
      call pgwindow(real(tmin),real(tmax),real(mmax),real(mmin))
      bb(1)=tmin
      bb(2)=tmax
      bb(3)=mmax
      bb(4)=mmin
C     range box, ticks and scale
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
C     create label for x-axis
      write(xlabel,500) "BJD-",MOSTtime+2400000.0
 500  format(A4,F16.8)
c      call pglabel(xlabel,"mmag","")
      call pgptxt(real((tmin+tmax)/2.0),real(mmax+0.31*(mmax-mmin)),
     .  0.0,0.5,xlabel)
      call pgptxt(real(tmin-0.05*(tmax-tmin)),real((mmax+mmin)/2),
     .  90.0,0.5,"mmag")
      call pgsch(2.0)
C     convert *8 to *4 and mag to mmag for plotting
      do 11 i=1,npt 
        px(i)=real(time(i))
        py(i)=real(mag(i)*1000.0)
        pz(i)=real(merr(i)*1000.0)
 11   continue
      call pgsch(0.8)
      call pgpt(npt,px,py,sym)
c      call pgsch(2.9)
c      call pgerrb(6,npt,px,py,pz,1.0)
      

      return
      end

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotph(npt,time,mag,merr,period,x1,x2,pan,sym,
     .  nmax,px,py,pz,phase,toff,nerrplot,bbp)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Plots the phased light curve
C     npt - number of data points
C     time - times of observations
C     mag  - the observations
C     merr - observational error
C     period - fixed period of the data
C     x1,x2 - phase bounds to plot
C     panx,pany - pgplot panel to plot light curve
C     sym - plotting symbol to use
C     toff - phase offset to center plots
C     nerrplot - 0: plot errors, 1: do not plot errors
C     nmax,px,py,pz,phase - workspace variables
C     bbp : bounds used for plot
      implicit none
      integer npt,pan,i,sym,nmax,nerrplot
      double precision time(npt),mag(npt),period,phase(nmax),mmax,mmin,
     .  ax1,ax2,stdev,std,sigcut,mean,merr(npt),x1,x2,bbp(4)
      real px(nmax),py(nmax),pz(nmax),toff

      sigcut=10.0
      mmin= 99.9d30
      mmax=-99.9d30

      mean=0.0d0
      do 5 i=1,npt
         mean=mean+mag(i)
 5    continue
      mean=mean/dble(npt)

      std=stdev(npt,mag,mean)
c      call pgpanl(panx,pany)
      call panel(pan)
      call phasept(npt,time,phase,period,toff)
      do 10 i=1,npt
        if((phase(i).ge.x1).and.(phase(i).le.x2))then
            if((mag(i).gt.mmax).and.
     .          (abs(mag(i)-mean).lt.sigcut*std))mmax=mag(i)
            if((mag(i).lt.mmin).and.
     .          (abs(mag(i)-mean).lt.sigcut*std))mmin=mag(i)
        endif
 10   continue

      ax2=1000.0*(mmin-0.10*(mmax-mmin))
      ax1=1000.0*(mmax+0.10*(mmax-mmin))
      bbp(1)=x1
      bbp(2)=x2
      bbp(3)=ax1
      bbp(4)=ax2
      call pgwindow(real(x1),real(x2),real(ax1),real(ax2))
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
      if(pan.le.2) then
        call pgptxt(real((x1+x2)/2.0),real(ax1+0.31*(ax1-ax2)),
     .      0.0,0.5,"Phase")
        call pgptxt(real(x1-0.05*(x2-x1)),real((ax2+ax1)/2),
     .      90.0,0.5,"mmag")
      elseif(pan.eq.3) then
        call pgptxt(real((x1+x2)/2.0),real(ax1+0.12*(ax1-ax2)),
     .      0.0,0.5,"Phase")
        call pgptxt(real(x1-0.1*(x2-x1)),real((ax2+ax1)/2),
     .      90.0,0.5,"mmag")
      endif
      do 11 i=1,npt
        px(i)=real(phase(i))
        py(i)=real(mag(i)*1000.0)
        pz(i)=real(merr(i)*1000.0)
 11   continue
      call pgpt(npt,px,py,sym)
      if(nerrplot.eq.0) call pgerrb(6,npt,px,py,pz,1.0)

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotrvph(nptv,vtime,vel,verr,period,x1,x2,pan,sym,
     .  nmax,px,py,pz,phase,toff,nerrplot,bbp)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Plots the phased light curve
C     npt - number of data points
C     time - times of observations
C     mag  - the observations
C     merr - observational error
C     period - fixed period of the data
C     x1,x2 - phase bounds to plot
C     panx,pany - pgplot panel to plot light curve
C     sym - plotting symbol to use
C     toff - phase offset to center plots
C     nerrplot - 0: plot errors, 1: do not plot errors
C     nmax,px,py,pz,phase - workspace variables
C     bbp : bounds used for plot
      implicit none
      integer nptv,pan,i,sym,nmax,nerrplot
      double precision vtime(nptv),vel(nptv),period,phase(nmax),mmax,
     .  mmin,ax1,ax2,stdev,std,sigcut,mean,verr(nptv),x1,x2,bbp(4)
      real px(nmax),py(nmax),pz(nmax),toff

      sigcut=10.0
      mmin= 99.9d30
      mmax=-99.9d30

      mean=0.0d0
      do 5 i=1,nptv
         mean=mean+vel(i)
 5    continue
      mean=mean/dble(nptv)

      std=stdev(nptv,vel,mean)
c      call pgpanl(panx,pany)
      call panel(pan)
      call phasept(nptv,vtime,phase,period,toff)
      do 10 i=1,nptv
        if((phase(i).ge.x1).and.(phase(i).le.x2))then
            if((vel(i).gt.mmax).and.
     .          (abs(vel(i)-mean).lt.sigcut*std))mmax=vel(i)
            if((vel(i).lt.mmin).and.
     .          (abs(vel(i)-mean).lt.sigcut*std))mmin=vel(i)
        endif
 10   continue

c      mmin=-50.0
c      mmax= 50.0

      ax1=(mmin-0.10*(mmax-mmin))
      ax2=(mmax+0.10*(mmax-mmin))
      bbp(1)=x1
      bbp(2)=x2
      bbp(3)=ax1
      bbp(4)=ax2
      call pgwindow(real(x1),real(x2),real(ax1),real(ax2))
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
      if(pan.le.2) then
        call pgptxt(real((x1+x2)/2.0),real(ax1+0.31*(ax1-ax2)),
     .      0.0,0.5,"Phase")
        call pgptxt(real(x1-0.1*(x2-x1)),real((ax2+ax1)/2),
     .      90.0,0.5,"Velocity (m/s)")
      elseif((pan.eq.3).or.(pan.eq.4)) then
        call pgptxt(real((x1+x2)/2.0),real(ax1+0.12*(ax1-ax2)),
     .      0.0,0.5,"Phase")
        call pgptxt(real(x1-0.2*(x2-x1)),real((ax2+ax1)/2),
     .      90.0,0.5,"Velocity (m/s)")
      endif
      do 11 i=1,nptv
        px(i)=real(phase(i))
        py(i)=real(vel(i))
        pz(i)=real(verr(i))
 11   continue
      call pgpt(nptv,px,py,sym)
      if(nerrplot.eq.0) call pgerrb(6,nptv,px,py,pz,1.0)

      return
      end
