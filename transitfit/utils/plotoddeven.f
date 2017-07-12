CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotoddeven(npt,time,mag,merr,period,x1,x2,pan,sym,
     .  nmax,px,py,pz,phase,toff,nerrplot,bbp,tdepth,mag2,merr2,bins,
     .  oddeven,phase2,doe,nsensor,flag)
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
      integer npt,pan,i,sym,nmax,nerrplot,npoints,npt2,bins,
     .  oddeven(npt),j,oe(2),k,nsensor,flag(npt)
      double precision time(npt),mag(npt),period,phase(nmax),mmax,mmin,
     .  ax1,ax2,stdev,std,sigcut,mean,merr(npt),x1,x2,bbp(4),width,
     .  tdepth,mag2(npt),merr2(npt),phase2(npt),x1t,x2t,shift,doe
      real px(nmax),py(nmax),pz(nmax),toff,pm(2),rbb(4)
      character*80 label

      oe(1)=0 !even transits
      oe(2)=1 !odd transits
      pm(1)=-1.0d0
      pm(2)= 1.0d0

      sigcut=3.0
      mmin= 99.9d30
      mmax=-99.9d30

      mean=0.0d0
      do 5 i=1,npt
         mean=mean+mag(i)
c         mag2(i)=mag(i)
c         merr2(i)=merr(i)
 5    continue
      mean=mean/dble(npt)

      std=stdev(npt,mag,mean)
c      call pgpanl(panx,pany)
      call panel(pan)
      call phaseptoddeven(npt,time,phase,period,toff,oddeven)
      npoints=0
      do 10 i=1,npt
        if((phase(i).ge.x1).and.(phase(i).le.x2))then
            if((mag(i).gt.mmax).and.
     .          (abs(mag(i)-mean).lt.sigcut*std))mmax=mag(i)
            if((mag(i).lt.mmin).and.
     .          (abs(mag(i)-mean).lt.sigcut*std))mmin=mag(i)
            if(abs(mag(i)-mean).lt.sigcut*std)npoints=npoints+1
        endif
 10   continue
c      write(0,*) "npoints:",npoints
      
c      write(0,*) "mmin",mmin
      if(1.0d0-mmin.lt.tdepth) mmin=1.0d0-tdepth-std
c      write(0,*) "mmin",mmin

      ax2=(mmax+0.10*(mmax-mmin))
      ax1=(mmin-0.10*(mmax-mmin))
      x1t=0.75-2.0*(0.75-x1)
      x2t=0.75+2.0*(x2-0.75)
      shift=(0.75d0-x1)
      bbp(1)=x1t
      bbp(2)=x2t
      bbp(3)=ax1
      bbp(4)=ax2
      call pgwindow(real(x1t),real(x2t),real(ax1),real(ax2))
c      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
c      if(pan.le.2) then
        call pgptxt(real((x1t+x2t)/2.0),real(ax1-0.11*(ax2-ax1)),
     .      0.0,0.5,"Phase (hours)")
        call pgptxt(real(x1t-0.2*(x2t-x1t)),real((ax2+ax1)/2),
     .      90.0,0.5,"Flux")
c      elseif(pan.eq.3) then
c        call pgptxt(real((x1t+x2t)/2.0),real(ax1-0.12*(ax2-ax1)),
c     .      0.0,0.5,"Phase (hours)")
c        call pgptxt(real(x1t-0.2*(x2t-x1t)),real((ax2+ax1)/2),
c     .      90.0,0.5,"Flux")
c      endif
      
      do 13 k=1,2
        npt2=0
        do 11 i=1,npt
            if((oddeven(i).eq.oe(k)).and.(phase(i).gt.x1).and.
     .        (phase(i).lt.x2))then
                npt2=npt2+1
                px(npt2)=real(phase(i)+pm(k)*shift)
                py(npt2)=real(mag(i))
                pz(npt2)=real(merr(i))
            endif
 11     continue
        if(npoints.lt.1000) then 
            sym=17
        else
            sym=1
        endif
        call pgpt(npt2,px,py,sym)
        if(nerrplot.eq.0) call pgerrb(6,npt2,px,py,pz,1.0)
 
        call pgsci(2)
        call pgpt1(0.75+real(pm(k)*shift),real(mmin),18)
c        write(0,*) "mmin:",mmin
        call pgsci(1)
      
        npt2=0
        do 6 i=1,npt
            if((oddeven(i).eq.oe(k)).and.(phase(i).gt.x1).and.
     .        (phase(i).lt.x2))then
                npt2=npt2+1
                phase2(npt2)=phase(i)
                mag2(npt2)=mag(i)
                merr2(npt2)=merr(i)
            endif
 6      continue
        call binp(npt2,phase2,mag2,merr2,bins,flag)
        j=0
        do 12 i=1,npt2
            if((phase2(i).gt.x1).and.(phase2(i).lt.x2).and.
     ,       (flag(i).eq.0))then
                j=j+1
                px(j)=real(phase2(i)+pm(k)*shift)
                py(j)=real(mag2(i))
                pz(j)=real(merr2(i))
            endif
 12     continue
        npt2=j
        call pgsci(4)
        if(npt2.gt.0)then
            call pgsch(1.5)
            call pgpt(npt2,px,py,17)
            call pgsch(0.8)
            if(x1.gt.0.5) call pgerrb(6,npt2,px,py,pz,1.0)
        endif
        call pgsci(1)
      
        if(k.eq.1) write(label,506) "Odd"
        if(k.eq.2) write(label,506) "Even"
 506    format(A4)
        call pgptxt(real(0.75+pm(k)*shift),
     .      real(bbp(4)-0.1*(bbp(4)-bbp(3))),0.0,0.5,label)
      
 13   continue
 
      px(1)=real(x1t)
      px(2)=real(x2t)
      py(1)=real(1.0d0-tdepth)
      py(2)=real(1.0d0-tdepth)
      call pgsch(3.0)
      call pgsci(2)
c      call pgsls(4)
      call pgline(2,px,py)
c      call pgsls(1)
      call pgsci(1)
      call pgsch(0.8)
      
      do 14 i=1,4
        rbb(i)=real(bbp(i))
 14   continue
      
      call pgsch(0.7)
      if(doe.gt.0.0d0)then
        write(label,505) "Depth-sig=",doe
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.25*(rbb(2)-rbb(1)),
     .      rbb(4)+0.04*(rbb(4)-rbb(3)),
     .      0.0,0.5,label)
 505    format(A10,F4.1)
      endif
      write(label,507) "TDepth=",tdepth*1.0d6," ppm"
 507  format(A7,F7.1,A4) 
      if(nsensor.lt.2) call pgptxt(rbb(1)+0.75*(rbb(2)-rbb(1)),
     .  rbb(4)+0.04*(rbb(4)-rbb(3)),
     .  0.0,0.5,label)
  
      call pgsch(0.8)
      
      if(x1.ge.0.5)then
        width=(x2t-x1t)*period*24.0d0/2.0d0
        call pgwindow(real(-width),real(width),real(ax1),real(ax2))
      else
        call pgwindow(real(x1t),real(x2t),real(ax1),real(ax2))
      endif
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
      
 
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine phaseptoddeven(npt,time,phase,period,toff,oddeven)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     npt - number of data points
C     time - times of observations
C     mag - the observations
C     phase - returned phase of observation
C     period - fixed period for data
      implicit none
      integer npt,oddeven(npt),tn
      double precision time(npt),phase(npt),period,toff,poff

      integer i
      double precision temp

      poff=toff*period

      do 10 i=1,npt
         temp=time(i)+poff
C        Get the phase
        tn=int(temp/period)
        oddeven(i)=mod(tn,2)
c        write(0,*) oddeven(i)
c        read(5,*)
        phase(i)=temp/period-int(temp/period)
C        apply optional phase offset to make plot pretty
c         phase(i)=phase(i)+toff
C        make sure phase is between 0 and 1
        if(phase(i).lt.0.0) phase(i)=phase(i)+1.0
        if(phase(i).gt.1.0) phase(i)=phase(i)-1.0
 10   continue
      return
      end
