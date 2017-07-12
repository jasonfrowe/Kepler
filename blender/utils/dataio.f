CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exportfit(nfit,sol,serr,Dpvary,err,doe,toff,titles)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,i
      double precision sol(nfit),serr(nfit,2),Dpvary(nfit),toff,
     .  err(nfit),doe
      include 'titles.f'
      
      open(unit=11,file="newfit.dat")
      
      do 10 i=1,nfit
        write(11,500) titles(i),sol(i),serr(i,1),serr(i,2),Dpvary(i),
     .      err(i)
 10   continue
      write(11,500) "OFF",toff,doe
 500  format(A3,5(1X,1PE17.10))
 
      close(11)
      
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
      subroutine readkeplc(nunit,nmax,npt,dtime,mag,merr,itime,
     .  Keplertime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nmax,npt,i
      double precision dtime(nmax),mag(nmax),merr(nmax),itime(nmax),
     .  Keplertime,sec2day,mintime,Pi
     
      Pi=acos(-1.d0)   !Pi
            
C     if time=0, then time=MOSTtime (that does not mean drink up)
c      Keplertime=54900.0
C     BJD Feb 5,2010 correction
      Keplertime=54900.0d0
c      Keplertime=0.5d0


C     number of seconds in a day
      sec2day=86400.0d0
      
      i=1
      
      mintime=99.9d30
  9   continue
 10   read(nunit,*,err=9,end=20) dtime(i),mag(i),merr(i)
cc  reset times form arbitrary start of 0.0 to HJD of field center
cc  add 53.038152 to MJD center 1st exp
c        dtime(i)=dtime(i)+53.038152
cc  add HJD 5/2/09 adjustment for field center
c        dtime(i)=dtime(i)+0.50020-Keplertime
c  add time dependent change for field center
c        dtime(i)=dtime(i)+4.1d-5*(dtime(i)-53.0)
        dtime(i)=dtime(i)-Keplertime
c        dtime(i)=dtime(i)+sin(2.0d0*pi*dtime(i)/372.5d0-
c     .      1.1208d0)*0.00280758d0
        dtime(i)=dtime(i)+0.5d0 !MJD half day offset

        mintime=min(mintime,dtime(i))
        mag(i)=-2.5*log10(mag(i)+1.0d0)
        merr(i)=0.00005
        itime(i)=1800.0/sec2day
        i=i+1
      goto 10
 20   continue
        
      npt=i-1
      write(0,*)   "-------------------------"
      write(0,500) "Observations read: ",npt
      write(0,*) "Mintime: ",mintime
      write(0,*)   "-------------------------"
 500  format(1X,A19,I6)
 
c      Keplertime=Keplertime+mintime  !correct time=0 time
c      do 30 i=1,npt
c         dtime(i)=dtime(i)-mintime
c 30   continue
 
      return
      end

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
      MOSTtime=2451545.0d0
c      MOSTtime=45.0d0


C     number of seconds in a day
      sec2day=86400.0d0
      
      i=1

      mintime=99.9d30 !find the minimum time of observations
   9  continue
  10  read(nunit,*,err=9,end=20) dtime(i),mag(i),merr(i),sky,xc,yc,fx,
     .      fy,fxy,tboard,mfield,ftotal,aflux,nap,itime(i),skystd
         itime(i)=30.0*60.0
         mintime=min(dtime(i),mintime)
         itime(i)=itime(i)/sec2day
c         if(mfield.lt.20000.0d0) goto 10 !cut SAA passages
         merr(i)=abs(merr(i))
         if(merr(i).le.0.0d0) goto 10 !make sure errors are positive
         rn=ntest(mag(i))
         if(rn.eq.1) goto 10  !test for NaNs
         rn=ntest(sky)
         if(rn.eq.1) goto 10  !test for NaNs
        i=i+1
      goto 10
 20   continue

      npt=i-1
      write(0,*)   "-------------------------"
      write(0,500) "Observations read: ",npt
      write(0,*) "Mintime: ",mintime
      write(0,*)   "-------------------------"
 500  format(1X,A19,I6)     
      
      MOSTtime=MOSTtime+mintime  !correct time=0 time
      do 30 i=1,npt
         dtime(i)=dtime(i)-mintime
 30   continue

 510  format(F13.8,1X,2(F9.6,1X),F8.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F9.2),1X,F8.3,1X,F6.2,1X,F8.4)

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