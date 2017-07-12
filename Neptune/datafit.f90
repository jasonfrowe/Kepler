program datafit
use precision
implicit none
integer :: nmax,i,filestatus,nunit,npt,nplot,np,nfit,nth,j,n,seed
integer, dimension(3) :: now
integer, allocatable, dimension(:) :: npix
real, allocatable, dimension(:) :: x,y
real(double) :: kmag,func,kernel
real(double), allocatable, dimension(:) :: time,flux,sky,xpos,ypos,xcoo,&
   ycoo,xp,yp,ep,sol,diff,th,xpos2,ypos2,mu
real(double), allocatable, dimension(:,:) :: Kr,Kri,Ryp
character(80) :: filename,fitsname

interface
   subroutine plotter(nplot,x,y)
      use precision
      implicit none
      integer, intent(inout) :: nplot
      real, dimension(:), intent(inout) :: x,y
   end subroutine plotter
end interface
interface
   subroutine fitdata(np,xp,yp,ep,xpos,ypos,nfit,sol)
      use precision
      implicit none
      integer, intent(inout) :: np,nfit
      real(double), dimension(:), intent(inout) :: xp,yp,ep,sol,xpos,ypos
   end subroutine fitdata
end interface
interface
   function inv(A) result(Ainv)
      use precision
      real(double), dimension(:,:), intent(in) :: A
      real(double), dimension(size(A,1),size(A,2)) :: Ainv
   end function inv
end interface

if(iargc().lt.1)then
   write(0,*) "Usage: datafit phot.dat"
   stop
endif

call getarg(1,filename)  !get filename for photometry

nunit=10 !unit number for file list
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

nmax=100000
allocate(time(nmax),flux(nmax),sky(nmax),xpos(nmax),ypos(nmax),         &
   npix(nmax),xcoo(nmax),ycoo(nmax))
i=1
do
   read(nunit,*,iostat=filestatus) time(i),flux(i),sky(i),kmag,xpos(i), &
      ypos(i),npix(i),fitsname,xcoo(i),ycoo(i)
   if(filestatus == 0) then
      i=i+1 !count number of files
      cycle
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!!",filename
      write(0,900) "iostat: ",filestatus,i
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit)
npt=i-1
write(0,*) "npt: ",npt

!call pgopen('?')
call pgopen('/xserve')
call pgask(.false.)
call pgpage()
call PGPAP ( 8.0 ,1.0)
call pgsch(1.2)
call pgslw(2)
call pgvport(0.15,0.85,0.15,0.85)


allocate(x(npt),y(npt),xp(npt),yp(npt),ep(npt))

np=0
do i=1,npt
   if((xcoo(i).gt.-0.1).and.(xcoo(i).lt.0.31).and.(ycoo(i).gt.-0.12)    &
    .and.(ycoo(i).lt.0.25))then
      np=np+1
      xp(np)=time(i)
      yp(np)=flux(i)
      ep(np)=sqrt(flux(i))
   endif
enddo

nplot=np
x=real(xp)
y=real(yp)
call plotter(nplot,x,y)

call position(np,xp,xpos,ypos)
call pixelfix(np,xp,yp,xpos,ypos)


nfit=4
allocate(sol(nfit))
sol(1)=8.05535e+06
sol(2)=-4722.27
sol(3)=62252.3
sol(4)=-66929.2

!Kernel guesses
nth=12
allocate(th(nth))
th=0.0
th(1)=0.0!0.5d6 !amplitude of the linear trend
th(2)=4.0!40.0d0 !scale (40 days)

th(3)=0.0!2.618400E+04 !amplitude of periodic trend
th(4)=1000.0 !decay time
th(5)=1.0 !smoothness
th(6)=1.0/0.1427106E+01 !period

th(7)=0.0!100.0 !amplitude of irregularities
th(8)=2.0 !scale days
th(9)=4.0 !shape parameter

th(10)=1000.0 !amplitude of correlated data
th(11)=1.0   !time scale
th(12)=sqrt(8.0d6) !photometric error

!allocate(Kr(np,np))
!write(0,*) "making kernel"
!do i=1,np
!   do j=1,np
!      Kr(i,j)=Kernel(nth,th,xp(i),xp(j))
!   enddo
!enddo
!write(0,*) "done.. "
!do i=1,10
!   write(0,501) (Kr(i,j),j=1,10)
!enddo
!allocate(Kri(np,np))
!write(0,*) "Inverting.."
!Kri=inv(Kr)
!write(0,*) "done.. "
!do i=1,10
!   write(0,501) (Kri(i,j),j=1,10)
!enddo

!allocate(mu(np))
!do i=1,np
!   mu(i)=func(nfit,sol,xpos(i))
!enddo

!call itime(now)
!seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
!n=1
!allocate(ryp(np,n))
!write(0,*) "make random set"
!call multinormal_sample(np,n,Kr,mu,seed,ryp)


!call fitdata(np,xp,yp,ep,xpos,ypos,nfit,sol)
!allocate(diff(np))
!do i=1,npt
!   diff(i)=flux(i)-func(nfit,sol,time(i))
!   write(6,500) time(i),diff(i),sqrt(flux(i)),xpos(i),ypos(i)
!enddo
500 format(108(1X,1PE17.10))
501 format(108(1X,1PE10.3))

!call pgpage()
!y=real(diff)
!call plotter(nplot,x,y)

do i=1,np
   write(6,500) xp(i),yp(i),xpos(i),ypos(i)
enddo

call pgclos()

end program datafit

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine position(np,time,xpos,ypos)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: np,i,j
integer, parameter :: nfit=5
real(double), dimension(nfit) :: xsol,ysol
real(double), dimension(np) :: time,xpos,ypos

xsol(1)=-239.532
xsol(2)= 17.7123
xsol(3)=-0.150237
xsol(4)=-0.00147985
xsol(5)=7.69349e-06

ysol(1)=122.219
ysol(2)=-3.06133
ysol(3)=-0.0267141
ysol(4)=0.00143165
ysol(5)=-8.22201e-06

do i=1,np
   xpos(i)=xsol(1)
   ypos(i)=ysol(1)
   do j=2,nfit
      xpos(i)=xpos(i)+xsol(j)*time(i)**(j-1)
      ypos(i)=ypos(i)+ysol(j)*time(i)**(j-1)
   enddo
enddo

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function kernel(nth,th,x1,x2)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nth,i
real(double) :: x1,x2,kernel,r2,r,Pi,delta
real(double), dimension(4) :: k
real(double), dimension(nth) :: th

Pi=acos(-1.d0)!define Pi and 2*Pi

r =(x1-x2)
r2=r*r
k(1)=th(1)*th(1)*exp(-r2/(2.0d0*th(2)*th(2)))
k(2)=th(3)*th(3)*exp(-r2/(2.0d0*th(4)*th(4))-th(5)*th(5)*sin(Pi*r/th(6))**2.0d0)
k(3)=th(7)*th(7)*(1.0d0+r2/(2.0d0*th(8)*th(9)*th(9)))**(-1.0d0*th(8))
k(4)=th(10)*th(10)*exp(-r2/(2.0d0*th(11)*th(11)))+th(12)*th(12)*delta(x1,x2)

kernel=k(1)+k(2)+k(3)+k(4)

!write(0,*) " "
!write(0,500) x1,x2,r
!write(0,500) (k(i),i=1,4),kernel
!500 format(108(1X,1PE17.10))
!read(5,*)

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function delta(x1,x2)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
real(double) :: delta,x1,x2

if(abs(x1-x2).lt.1.0d-5)then
   delta=1.0d0
else
   delta=0.0d0
endif

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fitdata(np,xp,yp,ep,xpos,ypos,nfit,sol)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittingmod
implicit none
integer :: lwa,info,n,m,nfit,i
integer, target :: np
integer, allocatable, dimension(:) :: iwa
real(double) :: tollm
real(double), dimension(:), target :: xp,yp,ep,sol,xpos,ypos
real(double), allocatable, dimension(:) :: wa,fvec,sol2
external fcn

write(0,*) "Fitting time"
n=nfit !number of variables to fit
m=np !number of data points
allocate(sol2(n))
sol2=sol

np2 => np
xp2 => xp
yp2 => yp
ep2 => ep
xpos2 => xpos
ypos2 => ypos


lwa=m*n+5*m*n
allocate(fvec(m),iwa(n),wa(lwa))
tollm=1.0d-12
!!call lmdif1(fcn,m,n,x,fvec,tollm,info,iwa,wa,lwa)
write(0,*) (sol2(i),i=1,n)
write(0,*) "Calling lmdif1"
!!call lmdif1(fcn,m,n,sol2,fvec,tollm,info,iwa,wa,lwa)
write(0,*) "info: ",info
write(0,*) (sol2(i),i=1,n)

sol=sol2

return
end subroutine fitdata

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn(npt,n,x,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittingmod
implicit none
integer npt,n,iflag,i
real(double) :: f,func
real(double), dimension(n) :: x
real(double), dimension(npt) :: fvec

!do i=1,npt
!   f=func(n,x,xp2(i))
!   fvec(i)=(yp2(i)-f)/ep2(i)
!enddo

return
end subroutine fcn

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function func(nfit,sol,xp)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nfit,i
real(double) :: xp,func,pi,tPi,X
real(double), dimension(nfit) :: sol
Pi=acos(-1.d0)!define Pi and 2*Pi
tPi=2.0d0*Pi

!func=sol(1)+sol(2)*xp
!func=func+sol(4)*cos(tPi*sol(3)*xp+sol(5))
!func=func+sol(7)*cos(tPi*sol(6)*xp+sol(8))

X=xp-floor(xp)

func=sol(1)
do i=2,nfit
   func=func+sol(i)*X**(i-1)
enddo

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine plotter(nplot,x,y)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nplot
real, allocatable, dimension(:) :: rj
real, dimension(:) :: x,y

allocate(rj(4))
rj(1)=minval(x(1:nplot))
rj(2)=maxval(x(1:nplot))
rj(3)=minval(y(1:nplot))
rj(4)=maxval(y(1:nplot))

call pgwindow(rj(1),rj(2),rj(3),rj(4))
call pgbox('BCNTS1',0.0,0,'BCNTS',0.0,0)
call pglabel("Time (days)","Counts (e-)","")

call pgpt(nplot,x,y,1)

return
end subroutine plotter
