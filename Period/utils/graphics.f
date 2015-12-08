CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotdata(npt,time,mag,merr,x1,x2,panx,pany,MOSTtime,
     .  title)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine will plot the data (ie. time vs mag)
C     integer npt: Total number of data points
C     real time(npt): The x-axis data points
C     real mag(npt): The y-axis data points
C     real merr(npt): Error associated with mag()
C     real x1,x2: range of data to plot
C     panx,pany: PGPLOT panel to plot data in

      implicit none
      integer npt,panx,pany,i,sym,nmax
      parameter(nmax=2000000)
      real px(nmax),py(nmax),bb(4),pz(nmax)
      double precision time(npt),mag(npt),merr(npt),x1,x2,mmin,mmax,
     .  tmin,tmax,MOSTtime
      character*80 xlabel,title

C     Plotting size
      sym=1
c      if(npt.gt.1000) sym=-1
     
      call pgpanl(panx,pany)

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

      if(x1.ne.-1.0) then
         tmin=x1
      endif
      if(x2.ne.-1.0) then
         tmax=x2
      endif

c      mmax=0.022
c      mmin=-0.008
c      mmax=22.0
c      mmin=-5.0
      write(6,*) "mmin, mmax",mmin,mmax
      bb(1)=real(tmin)
      bb(2)=real(tmax)
      bb(3)=real(mmin-0.1*(mmax-mmin))
      bb(4)=real(mmax+0.1*(mmax-mmin))
      call pgwindow(bb(1),bb(2),bb(3),bb(4))
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
      write(xlabel,500) "HJD-",MOSTtime
 500  format(A4,F16.8)
c      call pglabel(xlabel,"mmag","")
      call pgptxt((bb(1)+bb(2))/2.0,bb(3)-0.27*(bb(4)-bb(3))
     .  ,0.0,0.5,xlabel)
      call pgptxt((bb(1)+bb(2))/2.0,bb(4)+0.05*(bb(4)-bb(3))
     .  ,0.0,0.5,title)
      call pgptxt(bb(1)-0.05*(bb(2)-bb(1)),(bb(3)+bb(4))/2,
     .  90.0,0.5,"ppt")
      call pgsch(2.0)
      do 11 i=1,npt
        px(i)=real(time(i))
        py(i)=real(mag(i))
        pz(i)=real(merr(i))
c        write(6,*) px(i),py(i)
 11   continue
      call pgsch(3.0)
      call pgpt(npt,px,py,sym)
      call pgsch(2.9)
c      call pgerrb(6,npt,px,py,pz,1.0)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotph(npt,time,mag,period,panx,pany)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Plots the phased light curve
C
C     npt - number of data points
C     time - times of observations
C     mag  - the observations
C     period - fixed period of the data
C     panx,pany - pgplot panel to plot light curve
      implicit none
      integer nmax
      parameter(nmax=2000000)

      integer npt,panx,pany,i
      real px(nmax),py(nmax),bb(4),ax1,ax2
      double precision time(npt),mag(npt),period,phase(nmax),mmax,mmin,
     .  temp,stdev,std,sigcut,mean

      mean=0.
      do 4 i=1,npt
         mean=mean+mag(i)
 4    continue
      mean=mean/real(npt)

      sigcut=5.0
      mmax=-10.0e10
      mmin= 10.0e10

      std=stdev(npt,mag,mean)
      call pgpanl(panx,pany)
      call phasept(npt,time,mag,phase,period)
      do 10 i=1,npt
         if((mag(i).gt.mmax).and.
     .        (abs(mag(i)-mean).lt.sigcut*std))mmax=mag(i)
         if((mag(i).lt.mmin).and.
     .        (abs(mag(i)-mean).lt.sigcut*std))mmin=mag(i)
 10   continue

      ax2=real(mmin-0.10*(mmax-mmin))
      ax1=real(mmax+0.10*(mmax-mmin))
      call pgwindow(0.0,1.0,ax2,ax1)
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel("Phase","magnitude"," ")
      call pgptxt(0.5,ax2-0.27*(ax1-ax2),0.0,0.5,"Phase")
      call pgptxt(-0.05,(ax1+ax2)/2,90.0,0.5,"ppt")
      do 11 i=1,npt
        px(i)=real(phase(i))
        py(i)=real(mag(i))
 11   continue
c      if (npt.lt.1000) then
         call pgsch(3.0)
         call pgpt(npt,px,py,17)
         call pgsch(1.0)
c      else
c         call pgpt(npt,px,py,-1)
c      endif

      return
      end
