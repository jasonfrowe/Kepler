program koi3138srad
use precision
implicit none
integer :: iargc,nunit,nsolMCMCmax,nsolpars,filestatus,i,j,nptspars,    &
 ntpars,npttpars,seed,k,nbin,nx,ny,ncol
integer, dimension(3) :: now
integer, allocatable, dimension(:,:) :: ia
real :: x,y,rbin,Red,Green,Blue,val,z1,z2,dx,dy,px,py
real, allocatable, dimension(:) :: rj,contours,tr
real, allocatable, dimension(:,:) :: array
real(double) :: dumd,ran2,Psec,M1,R1,aConst,asemi,Pi,Teff,srad,Rp,rdr
real(double), allocatable, dimension(:) :: tparssol
real(double), allocatable, dimension(:,:) :: sparssol
character(80) :: sMCMCfile,tMCMCfile

if(iargc().lt.2)then
   write(0,*) "Usage: koi3138srad transitmcmc stellarmcmc"
   write(0,*) " transitmcmc - MCMC analysis of transit"
   write(0,*) " stellarmcmc - MCMC analysis of stellar parameters"
   stop
endif

Pi=acos(-1.d0)   !Pi
aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0) !for calculating semi-major axis
call itime(now)
seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
dumd=ran2(-seed)

call getarg(2,sMCMCfile)  !get filename for stellarMCMC

nsolMCMCmax=100001 !initital guess at the number of chains
nsolpars=7 !number of parameters
allocate(sparssol(nsolMCMCmax,nsolpars))
nunit=10 !unit number for data spectrum
open(unit=nunit,file=sMCMCfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",sMCMCfile
   stop
endif

i=1
do
   if(i.gt.nsolMCMCmax)then !check that we have enough allocated memory
      write(0,*) "Critical Error: Increase nsolMCMCmax"
      stop
   endif
   read(nunit,*,iostat=filestatus) (sparssol(i,j),j=1,nsolpars)
   if(filestatus == 0) then
      i=i+1
      cycle
   elseif(filestatus == -1) then
      exit
   else
      write(0,*) "File Error!!"
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit)
nptspars=i-1  !number of chains read from stellar MCMC work
write(0,*) "nptspars: ",nptspars

!now we set up our plot.

call pgopen('?') !open plotting device
call pgpage()
call pgask(.false.)
call PGPAP ( 8.0 ,1.0) !use a square 8" across
!call pgsubp(1,4)
call pgpage()
call pgsch(1.5) !make the font a bit bigger
call pgslw(3)  !make the lines a bit thicker
call pgvport(0.15,0.85,0.15,0.85) !make room around the edges for labels

allocate(rj(4))
rj(1)=log10(0.2)! 0.0 !scale for plot
rj(2)=log10(5.0)!5.0
rj(3)=log10(0.4)!0.0
rj(4)=log10(2.0)!2.0
call pgwindow(rj(1),rj(2),rj(3),rj(4)) !plot scale
call pgbox("BCLNTS1",0.0,0,"BCLNTS1",0.0,0)
call pglabel("Incident Flux (S\d\(2284)\u)","Radius (R\d\(2284)\u)","")


nbin=100 !binned for 'hess-diagram'
rbin=real(nbin)
allocate(array(nbin,nbin))
array=0

ntpars=18
allocate(tparssol(ntpars))
call getarg(1,tMCMCfile)  !get filename for stellarMCMC

open(unit=nunit,file=tMCMCfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",tMCMCfile
   stop
endif
!call pgsci(15)
!call pgbbuf()
i=1
do
   read(nunit,*,iostat=filestatus) dumd,dumd,dumd,(tparssol(j),j=1,ntpars)
   k=ran2(seed)*(nptspars-1)+1 !random number to pick a stellar MC
   Psec=tparssol(10)*day !period in seconds
   M1=sparssol(k,1)*Msun !mass in KG
   R1=sparssol(k,4)*Rsun !radius in m
   Teff=sparssol(k,6) !Teff in K
   asemi=(M1)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst
   srad=sparssol(k,4)*sparssol(k,4)*(Teff/Tsun)**4.0*(AU/asemi)**2.0
   rdr=tparssol(12)
   Rp=sparssol(k,4)*rdr*109.17

   x=real(log10(srad))
   y=real(log10(Rp))
!   call pgpt1(x,y,-1)

   nx=int( (x-rj(1))/(rj(2)-rj(1))*rbin )
   ny=int( (y-rj(3))/(rj(4)-rj(3))*rbin )
!   write(0,*) x,y,nx,ny
   if( (nx.gt.0).and.(nx.le.nbin).and.(ny.gt.0).and.(ny.le.nbin) )then
      array(nx,ny)=array(nx,ny)+1.0
   endif

   if(filestatus == 0) then
      i=i+1
      cycle
   elseif(filestatus == -1) then
      exit
   else
      write(0,*) "File Error!!"
      write(0,900) "iostat: ",filestatus
      stop
   endif
enddo
close(nunit)
!call pgebuf()
npttpars=i-1  !number of chains read from stellar MCMC work
write(0,*) "npttpars: ",npttpars
!call pgsci(1)

ncol=64 !set up colour-map
do i=1,ncol
!   Red = REAL(I-1)/REAL(NCOL-1)*0.8 + 0.2
!   Green = MAX(0.0, 2.0*REAL(I-1-NCOL/2)/REAL(NCOL-1))
!   Blue = 0.2 + 0.4*REAL(NCOL-I)/REAL(NCOL)
   Red = 1.0 - real(i)/real(ncol)
   Green = 1.0 - real(i)/real(ncol)*0.5
   Blue = 1.0 - real(i)/real(ncol)*0.3
   CALL PGSCR(I+15, Red, Green, Blue)
enddo

do i=1,nbin
   do j=1,nbin
!      write(0,*) rj(3)+real(i)*(rj(4)-rj(3))/rbin
      px=10.0**(rj(1)+real(i+1)*(rj(2)-rj(1))/rbin)-                      &
         10.0**(rj(1)+real(i)*(rj(2)-rj(1))/rbin)
      py=10.0**(rj(3)+real(j+1)*(rj(4)-rj(3))/rbin)-                      &
         10.0**(rj(3)+real(j)*(rj(4)-rj(3))/rbin)
!      write(0,*) i,j,px,py
      array(i,j)=array(i,j)/px/py
!      read(5,*)
   enddo
enddo

z1=minval(array)
z2=maxval(array)
allocate(ia(nbin,nbin))
do i=1,nbin
   do j=1,nbin
      val=array(i,j)
      IA(i,j)=int((val-z1)/(z2-z1)*(NCOL-1))+16
      if(val.lt.z1) ia(i,j)=16
      if(val.gt.z2) ia(i,j)=ncol+15
   enddo
enddo

call pgpixl(ia,nbin,nbin,1,nbin,1,nbin,rj(1),rj(2),rj(3),rj(4))

allocate(contours(3),tr(6))
!do i=1,nbin
!   do j=1,nbin
!      rarray(i,j)=real(array(i,j))
!   enddo
!enddo
contours(1)=z1+(z2-z1)*0.6
contours(2)=z1+(z2-z1)*0.15
contours(3)=z1+(z2-z1)*0.02
tr(1)=rj(1)
tr(2)=(rj(2)-rj(1))/rbin
tr(3)=0.0
tr(4)=rj(3)
tr(5)=0.0
tr(6)=(rj(4)-rj(3))/rbin
call pgsci(14)
call pgcont(array,nbin,nbin,1,nbin,1,nbin,contours,-3,TR)
call pgsci(1)

call pgbox("BCLNTS1",0.0,0,"BCLNTS1",0.0,0)



dx=0.03
dy=-0.01
px=log10(1.0)
py=log10(1.0)
call pgsch(1.05)
call pgpt1(px,py,17) !Earth
call pgptxt(px+dx,py+dy,0.0,0.0,"Earth")

px=log10(0.4307)
py=log10(0.5317)
call pgpt1(px,py,17) !Mars
call pgptxt(px+dx,py+dy,0.0,0.0,"Mars")

px=log10(1.9113)
py=log10(0.95044)
call pgpt1(px,py,17) !Venus
call pgptxt(px+dx,py+dy,0.0,0.0,"Venus")

dx=0.04
dy=-0.025

px=0.320
py=1.11
call pgpt1(log10(px),log10(py),17) !K186
call pgerr1(1,log10(px),log10(py),log10((px+0.059)/px),0.0)
call pgerr1(3,log10(px),log10(py),log10((px+0.039)/px),0.0)
call pgerr1(2,log10(px),log10(py),log10((py+0.14)/py),0.0)
call pgerr1(4,log10(px),log10(py),log10((py+0.13)/py),0.0)
call pgptxt(log10(px)+dx,log10(py)+dy,0.0,0.0,"K186f")

px=1.59
py=1.74
call pgpt1(log10(px),log10(py),17) !K69c
call pgerr1(1,log10(px),log10(py),log10((px+0.88)/px),0.0)
call pgerr1(3,log10(px),log10(py),log10((px+0.42)/px),0.0)
call pgerr1(2,log10(px),log10(py),log10((py+0.34)/py),0.0)
call pgerr1(4,log10(px),log10(py),log10((py+0.20)/py),0.0)
call pgptxt(log10(px)+dx,log10(py)+dy,0.0,0.0,"K69c")

px=1.17
py=1.73
call pgpt1(log10(px),log10(py),17) !K62e
call pgerr1(1,log10(px),log10(py),log10((px+0.21)/px),0.0)
call pgerr1(3,log10(px),log10(py),log10((px+0.21)/px),0.0)
call pgerr1(2,log10(px),log10(py),log10((py+0.08)/py),0.0)
call pgerr1(4,log10(px),log10(py),log10((py+0.09)/py),0.0)
call pgptxt(log10(px)-dx,log10(py)+dy,0.0,1.0,"K62e")


px=0.41
py=1.42
call pgpt1(log10(px),log10(py),17) !K62f
call pgerr1(1,log10(px),log10(py),log10((px+0.08)/px),0.0)
call pgerr1(3,log10(px),log10(py),log10((px+0.07)/px),0.0)
call pgerr1(2,log10(px),log10(py),log10((py+0.06)/py),0.0)
call pgerr1(4,log10(px),log10(py),log10((py+0.08)/py),0.0)
call pgptxt(log10(px)-dx,log10(py)+dy,0.0,1.0,"K62f")


999 call pgclos()

end program koi3138srad
