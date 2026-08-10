Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readdata(nunit,nmax,npt,dtime,mag,merr,itime,MOSTtime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in data
C     npt - number of data points (returned)
C     time,mag,merr - the data and associated error.
C     itime - integration time (day) [same units as time]
      implicit none
      integer nmax,nunit
      integer npt,i,ntest,rn
      double precision dtime(nmax),mag(nmax),merr(nmax),sky,xc,yc,fx,fy,
     .  fxy,tboard,mfield,ftotal,aflux,nap,itime(nmax),skystd,
     .  mintime,MOSTtime,sec2day

C     if time=0, then time=MOSTtime (that does not mean drink up)
c      MOSTtime=2451545.0d0
      MOSTtime=51545.0d0


C     number of seconds in a day
      sec2day=86400.0d0
      
      i=1

      mintime=99.9d30 !find the minimum time of observations
   9  continue
  10  read(nunit,*,err=9,end=20) dtime(i),mag(i),merr(i)!,sky,xc,yc,fx,
C     .      fy,fxy,tboard,mfield,ftotal,aflux,nap,itime(i),skystd
         itime(i)=10.26
         mintime=min(dtime(i),mintime)
         itime(i)=itime(i)/sec2day
c         if(mfield.lt.20000.0d0) goto 10 !cut SAA passages
         mag(i)=10.0**(mag(i)/(-2.5d0))
         merr(i)=abs(merr(i))
         if(merr(i).le.0.0d0) goto 10 !make sure errors are positive
c         rn=ntest(mag(i))
c         if(rn.eq.1) goto 10  !test for NaNs
c         rn=ntest(sky)
c         if(rn.eq.1) goto 10  !test for NaNs
        i=i+1
      goto 10
 20   continue

      npt=i-1
      write(0,*)   "-------------------------"
      write(0,500) "Observations read: ",npt
      write(0,*) "Mintime: ",mintime
      write(0,*)   "-------------------------"
 500  format(1X,A19,I6)     
      
c      MOSTtime=MOSTtime+mintime  !correct time=0 time
c      do 30 i=1,npt
c         dtime(i)=dtime(i)-mintime
c 30   continue

 510  format(F13.8,1X,2(F9.6,1X),F8.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F9.2),1X,F8.3,1X,F6.2,1X,F8.4)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exportfit(nfit,nplanet,sol,serr,err,titles)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,i,nplanet
      double precision sol(nfit),serr(nfit,2),err(nfit)
      character*3 tit
      include 'titles.f'
      
      open(unit=11,file="newfit.dat")

      sol(1)=abs(sol(1))

      do 13 i=1,nplanet
         sol(11+10*(i-1))=abs(sol(11+10*(i-1)))
 13   continue

      do 10 i=1,nplanet*10+8
        if(i.le.8)then
            tit=titles(i)
        else
            tit=titles(i-10*((i-9)/10))
!            write(0,*) titles(i-10*((i-9)/10))
            write(tit(3:3),501) ((i-9)/10)+1
 501        format(I1)
!            write(0,*) i,i-10*((i-9)/10),tit
        endif
        write(11,500) tit,sol(i),serr(i,1),serr(i,2),err(i)
 10   continue
 500  format(A3,5(1X,1PE17.10))
 
      close(11)
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readkeplc(nunit,nmax,npt,dtime,flux,ferr,itime,
     .  Keplertime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nmax,npt,i
      double precision dtime(nmax),flux(nmax),ferr(nmax),itime(nmax),
     .  Keplertime,sec2day,mintime,it
      character*400 line
      character*20 zeros


      write(zeros, '(*(I2))') (0, i=1, 4) !write a bunch of zeros
            
C     if time=0, then time=MOSTtime (that does not mean drink up)
      Keplertime=54900.0
c      Keplertime=0.5d0

C     number of seconds in a day
      sec2day=86400.0d0
      
      i=1
      
      it=0.0d0 !default, incase column is missing.

      mintime=99.9d30
  9   continue
 10   read(nunit,'(A)',end=20) line
      line=trim(line) // zeros  !dump zeros to the end.. 
      read(line,*,err=9,end=20) dtime(i),flux(i),ferr(i),it
        dtime(i)=dtime(i)-Keplertime
        dtime(i)=dtime(i)+0.5d0 !MJD half day offset
        mintime=min(mintime,dtime(i))
        flux(i)=flux(i)+1.0!-2.5*log10(mag(i)+1.0d0)
c        ferr(i)=0.00005
        if ((it.lt.1.0d-7).and.(it.gt.-1.0d-7)) then
          !itime(i)=1765.5/sec2day
          itime(i)=1764.944/sec2day
        elseif (it.lt.0.0) then
          itime(i)=58.85/sec2day !short cadence
        else
          itime(i)=it
        endif
c        itime(i)=0.1d0/sec2day
c         itime(i)=itime(i)/sec2day
        !write(0,*) dtime(i),flux(i),ferr(i),itime(i),it
        !read(5,*)
        i=i+1
      goto 10
 20   continue
        
      npt=i-1
      write(0,*)   "-------------------------"
      write(0,500) "Observations read: ",npt
c      write(0,*) "Mintime: ",mintime
      write(0,*)   "-------------------------"
 500  format(1X,A19,I8)
 
c      Keplertime=Keplertime+mintime  !correct time=0 time
c      do 30 i=1,npt
c         dtime(i)=dtime(i)-mintime
c 30   continue
 
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readrv(nunit,nmax,nptv,vtime,vel,verr,vetime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nptv,i,nmax
      double precision vtime(nmax),vel(nmax),verr(nmax),vetime(nmax),
     .  RVtime,sec2day
      character dumc

C     number of seconds in a day
      sec2day=86400.0d0
C     RV offset
      RVtime=4900.0d0
      
C     First 5 lines are comments
      do 10 i=1,5
        read(nunit,*) dumc
 10   continue
 
      i=1
 11   read(nunit,*,end=12) vtime(i),vel(i),verr(i) 
        vtime(i)=vtime(i)-RVtime
        vetime(i)=1800.0/sec2day
        i=i+1
      goto 11
 12   continue
      nptv=i-1
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function ntest(rr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     this routine is used to test for NANs and other nasty numbers.
      implicit none
      double precision rr,thi,tlow

C     this defines the liberal numerical boundary of anything that 
C     might be exciting in MOST data

      thi=  99.9d30
      tlow=-99.9d30

      if((rr.gt.tlow).and.(rr.lt.thi))then
         ntest=0
      else
         ntest=1
      endif

      return
      end      
