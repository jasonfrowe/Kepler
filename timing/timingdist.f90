program timingdist
implicit none
integer nmax,nunit,ndist,filestatus,i,nkid,j,seed,neccn
integer, dimension(3) :: now
integer, allocatable, dimension(:) :: dkid,kid
real :: a,ec,w,theta,timing,ran2,dumr,mean
real, allocatable, dimension(:) :: din,dist,period,dt,eccn,eprb
real, parameter :: pi = 3.1415926535897932
character(80) :: dfile,kidfile,eccnfile

dfile="KQ17distances.dat"  !file with distances
kidfile="Q1Q17KOIKID.dat"  !file with KIDs of planet candidates
eccnfile="eccndisp.dat"    !file with eccentricity probability dist.

nmax=200000
allocate(dkid(nmax),din(nmax))
!din - distances that are read in

nunit=10
open(unit=nunit,file=dfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",dfile
   stop
endif

i=1
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) dkid(i),din(i)
   if(filestatus == 0) then
      i=i+1
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
ndist=i-1
write(0,*) "ndist: ",ndist

allocate(kid(nmax),period(nmax))
!kid contains KIDs for planet candidates
!period contains orbital period for planet candidates

nunit=10
open(unit=nunit,file=kidfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",kidfile
   stop
endif

i=1
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) kid(i),period(i)
   if(filestatus == 0) then
      i=i+1
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      stop
   endif
enddo
close(nunit) !close file
nkid=i-1
write(0,*) "nkid: ",nkid

!read in eccentricity distribution
nunit=10
open(unit=nunit,file=eccnfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",kidfile
   stop
endif

neccn=100+1 !number of entries in eccnfile, with extra buffer to avoid errors
allocate(eccn(neccn),eprb(neccn)) !allocate space

i=1
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) eccn(i),eprb(i)
   if(filestatus == 0) then
      i=i+1
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      stop
   endif
enddo
close(nunit) !close file
neccn=i-1
write(0,*) "neccn: ",neccn


allocate(dist(nmax))
!dist contains dist for each KID
dist=0.0 !allocate to zero

do i=1,nkid
   do j=1,ndist
      if(kid(i).eq.dkid(j))then
         dist(i)=din(j)
      endif
   enddo
!   write(0,*) kid(i),dist(i),period(i)
enddo
write(0,*) "Distances allocated"

!set up seed for random number generator
call itime(now)
seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
dumr=ran2(-seed)

allocate(dt(nmax))
a=1.0 !AU
ec=0.0 !eccentricity
w=0.0 !omega
theta=0.0 !allignment
mean=0.0
do i=1,nkid
   if(dist(i).gt.0.0)then
      dt(i)=timing(a,dist(i),ec,w,theta,period(i))
   else
      dt(i)=0.0d0
   endif
   mean=mean+dt(i)
   dt(i)=log10(dt(i)+1.0e-8)
!   write(6,*) kid(i),dist(i),period(i),dt(i)
enddo
mean=mean/real(nkid)
write(0,*) "Mean: ",mean

!set up PGPLOT
!open PGPLOT device
call pgopen('?')  !'?' lets the user choose the device.
call PGPAP (8.0 ,1.0) !use a square 8" across
!call pgsubp(2,2)
!call pgpage() !create a fresh page
call pgslw(3) !thicker lines

!call pgpanl(2,2)
call pghist(nkid,dt,-6.0,0.05,50,0)
call pglabel("log dt (sec)","# of Planetary Candidates","")
!
!call pgpanl(1,1)
!! random allignment
!do i=1,nkid
!   theta=0.0!2.0*pi*ran2(seed)
!   ec=ran2(seed)
!   w=2.0*pi*ran2(seed)
!   if(dist(i).gt.0.0)then
!      dt(i)=timing(a,dist(i),ec,w,theta,period(i))
!   else
!      dt(i)=0.0d0
!   endif
!   dt(i)=log10(dt(i)+1.0e-8)
!!   write(6,*) kid(i),dist(i),period(i),dt(i)
!enddo

!call pghist(nkid,dt,-6.0,0.05,50,0)
!call pglabel("log dt (sec)","# of PCs","")

call pgclos()


end program timingdist

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
real function timing(a,d,e,w,th,P)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
implicit none
real ::  a,d,e,w,th,P
real :: AU,PC,PS
real, parameter :: pi = 3.1415926535897932

AU = 1.496e11 !1 AU (m)
PC = 3.086e16 !1 pc (m)
PS = 86400.0  !number of seconds in a day

timing = a*AU*P*PS*sqrt(1.0-e*e)/(2.0*pi*d*PC*(1+e*sin(w))*cos(th))

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION ran2(idum)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,    &
       IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,               &
       IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END

