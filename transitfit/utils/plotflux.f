CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotflux(npt,time,mag,merr,x1,x2,pan,MOSTtime,nmax,px,
     .  py,pz,bb,nfit,sol,id,id2,nsensor,kID,Kmag,Teff,logg,rad,tdepth)
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
      integer npt,i,sym,nmax,pan,nfit,id,id2,nsensor,kID,Teff
      real px(nmax),py(nmax),pz(nmax)
      double precision time(npt),mag(npt),merr(npt),x1,x2,mmin,mmax,
     .  tmin,tmax,bb(4),tmp(2),mean,stdev,std,sigcut,sol(nfit),T0,lpt,
     .  Kmag,logg,rad,tdepth
      double precision MOSTtime
      character*80 xlabel,label

C     If we have few observations, use big dots, otherwise, small dots
      sym=17
      if(npt.gt.100) sym=-1

C     select panel for plot
c      call pgpanl(panx,pany)
      call panel(pan)
      call pgsch(0.8)

      sigcut=3.0

      mean=0.0d0
      do 5 i=1,npt
         mean=mean+mag(i)
 5    continue
      mean=mean/dble(npt)

      std=stdev(npt,mag,mean)
      
C     find the boundaries of the data
      mmin=mag(1)
      mmax=mag(1)
      tmin=time(1)
      tmax=time(1)
      do 10 i=2,npt
         if((mag(i).gt.mmax).and.
     .      (abs(mag(i)-mean).lt.sigcut*std)) mmax=mag(i)
         if((mag(i).lt.mmin).and.
     .      (abs(mag(i)-mean).lt.sigcut*std)) mmin=mag(i)
         if(time(i).gt.tmax) tmax=time(i)
         if(time(i).lt.tmin) tmin=time(i)
 10   continue
 
c      write(0,*) "mmin",mmin
      if(1.0d0-mmin.lt.tdepth) mmin=1.0d0-tdepth-std
c      write(0,*) "mmin",mmin

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
      lpt=mmin
      tmp(1)=mmin-0.1*(mmax-mmin)
      tmp(2)=mmax+0.1*(mmax-mmin)
      mmin=tmp(1)
      mmax=tmp(2)
      call pgwindow(real(tmin),real(tmax),real(mmin),real(mmax))
      bb(1)=tmin
      bb(2)=tmax
      bb(3)=mmin
      bb(4)=mmax
C     range box, ticks and scale
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
c1      call pgbox('BCTS1',0.0,0,'BCNTSV1',0.0,0)
C     create label for x-axis
      write(xlabel,510) "BJD-",MOSTtime+2400000.0
 510  format(A4,F9.1)
c      call pglabel(xlabel,"mmag","")
      call pgptxt(real((tmin+tmax)/2.0),real(mmin-0.31*(mmax-mmin)),
     .  0.0,0.5,xlabel)
      call pgptxt(real(tmin-0.1*(tmax-tmin)),real((mmax+mmin)/2),
     .  90.0,0.5,"Flux")
      call pgsch(2.0)
C     convert *8 to *4 and mag to mmag for plotting
      do 11 i=1,npt 
        px(i)=real(time(i))
        py(i)=real(mag(i))
        pz(i)=real(merr(i))
 11   continue
      call pgsch(0.8)
      call pgpt(npt,px,py,sym)
c      call pgsch(2.9)
c      call pgerrb(6,npt,px,py,pz,1.0)

      call pgsci(2)
      T0=sol(7)-int((sol(7)-Tmin)/sol(5))*sol(5)
      do 12 i=1,int((Tmax-T0)/sol(5))+1
        px(1)=real(T0)+real(i-1)*real(sol(5))
        py(1)=lpt
        call pgpt(1,px,py,13)
 12   continue
      
      call pgsci(1)
      
      call pgsch(0.7)
      
      write(label,500) "Kepler:",id,".0",id2
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
      
      
      call pgsch(0.8)

      return
      end
