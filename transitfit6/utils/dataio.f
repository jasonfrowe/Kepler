CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readkeplc(nunit,nmax,npt,dtime,flux,ferr,itime,
     .  Keplertime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nmax,npt,i
      double precision dtime(nmax),flux(nmax),ferr(nmax),itime(nmax),
     .  Keplertime,sec2day,mintime
            
C     if time=0, then time=MOSTtime (that does not mean drink up)
      Keplertime=54900.0
c      Keplertime=0.5d0

C     number of seconds in a day
      sec2day=86400.0d0
      
      i=1
      
      mintime=99.9d30
  9   continue
 10   read(nunit,*,err=9,end=20) dtime(i),flux(i),ferr(i)
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
        flux(i)=flux(i)+1.0!-2.5*log10(mag(i)+1.0d0)
c        ferr(i)=0.00005
        itime(i)=1765.5/sec2day
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
      subroutine exportfit(nfit,nplanet,sol,serr,err,titles)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,i,nplanet
      double precision sol(nfit),serr(nfit,2),err(nfit)
      character*3 tit
      include 'titles.f'
      
      open(unit=11,file="newfit.dat")

      sol(1)=abs(sol(1))      
      do 10 i=1,nplanet*11+15
        if(i.le.15)then
            tit=titles(i)
        else
            tit=titles(i-11*((i-15-1)/11))
            write(tit(3:3),501) ((i-15-1)/11)+1
 501        format(I1)
c            write(0,*) i,i-10*((i-9)/10)
        endif
        write(11,500) tit,sol(i),serr(i,1),serr(i,2),err(i)
 10   continue
 500  format(A3,5(1X,1PE17.10))
 
      close(11)
      
      return
      end