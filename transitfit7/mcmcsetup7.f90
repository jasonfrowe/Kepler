program mcmcsetup7
use precision
use fittingmod
implicit none
integer nmax,i,j,nplanet,ii
integer, target :: nbodies,npt
real(double) :: drerr,dfridr,h,mdiff,Pi,LU2,rstar,length,chi1,chi2,dchi
real(double), target :: tol
real(double), allocatable, dimension(:) :: soldiff,ydiff,y3,m3,sol3,ans
real(double), allocatable, dimension(:), target :: time,flux,ferr,itime,sol,y, &
   err,yerr
real(double), allocatable, dimension(:,:), target :: serr,yserr
character(80) :: photfile,inputsol
external func

interface
   subroutine readkeplc(photfile,npt,time,flux,ferr,itime)
      use precision
      implicit none
      character(80), intent(in) :: photfile
      integer, intent(out) :: npt
      real(double), dimension(:), intent(inout) :: time,flux,ferr,itime
   end subroutine readkeplc
end interface

interface
   subroutine calcnbodies(inputsol,nbodies)
      use precision
      implicit none
      character(80), intent(in) :: inputsol
      integer, intent(out) :: nbodies
   end subroutine calcnbodies
end interface

interface
   subroutine readinputsol(inputsol,nbodies,y,yserr,yerr,sol,serr,err)
      use precision
      implicit none
      character(80), intent(in) :: inputsol
      integer, intent(in) :: nbodies
      real(double), dimension(:), intent(out) :: y,yerr,sol,err
      real(double), dimension(:,:), intent(out) :: yserr,serr
   end subroutine readinputsol
end interface

interface
   subroutine exportfit(nplanet,sol,serr,err,y,yerr,yserr)
      use precision, only: double
      implicit none
      integer, intent(in) :: nplanet
      real(double), dimension(:), intent(in) :: sol,err,y,yerr
      real(double), dimension(:,:), intent(in) :: serr,yserr
   end subroutine exportfit
end interface

interface
   subroutine lcmodel(nbodies,npt,tol,y,sol,time,itime,ans)
      use precision
      implicit none
      integer, intent(in) :: nbodies,npt
      real(double), intent(in) :: tol
      real(double), dimension(:), intent(inout) :: y,ans
      real(double), dimension(:), intent(in) :: sol,time,itime
   end subroutine lcmodel
end interface

Pi=acos(-1.d0)!define Pi
LU2=LU*LU

tol=1.0d-12  !default tolerance level

if(iargc().lt.2) then
   write(0,*) "Usage: mcmcsetup7 <photometry> <inputsol>"
   stop
endif

call getarg(1,photfile)
nmax=100000
allocate(time(nmax),flux(nmax),ferr(nmax),itime(nmax))
call readkeplc(photfile,npt,time,flux,ferr,itime)
write(0,*) "Number of observations read: ",npt

call getarg(2,inputsol)
call calcnbodies(inputsol,nbodies)
allocate(m(nbodies),mserr(nbodies,2),merr(nbodies),y(nbodies*6), &
   yserr(nbodies*6,2),yerr(nbodies*6),sol(8+nbodies), &
   serr(8+nbodies,2),err(8+nbodies))
write(0,*) "nbodies: ",nbodies
call readinputsol(inputsol,nbodies,y,yserr,yerr,sol,serr,err)
do i=1,nbodies
   write(0,500) m(i),(y(j),j=6*i-5,6*i)
enddo
write(0,*) " "
500  format(28(1PE10.3,1X))

!assign pointers from module that is shared with FCN
tol2 => tol
nbodies2 => nbodies
sol2 => sol
y2 => y
serr2 => serr
yserr2 => yserr
time2 => time
flux2 => flux
ferr2 => ferr
itime2 => itime
npt2 => npt

allocate(y3(nbodies2*6),sol3(8+nbodies2),m3(nbodies2),ans(npt))

m3=m
y3=y
sol3=sol
call lcmodel(nbodies,npt,tol,y3,sol3,time,itime,ans)
m=m3 !store mass values to module
chi1=0
do i=1,n
   chi1=chi1+(flux(i)-ans(i))*(flux(i)-ans(i))/(ferr(i)*ferr(i))
enddo

isol=0 !tells func which variable is being tweeked
im=0
iy=0
allocate(soldiff(9))
soldiff(1)=min(0.1,0.1*sol(1)) !mean stellar density
soldiff(2)=0.001 !limb
soldiff(3)=0.001
soldiff(4)=0.001
soldiff(5)=0.001
soldiff(6)=0.001
soldiff(7)=0.5
soldiff(8)=1.0d-6
do i=1,8
   isol=i
   if(serr(i,2).ne.0.0d0)then
      h=soldiff(i)
      serr(i,2)=1.0d0/abs(dfridr(func,sol(i),h,drerr))
   endif
   write(0,500) sol(i),serr(i,2)

   dchi=0.0d0
   do while((dchi.lt.0.5).or.(dchi.lt.2.0))
      m3=m
      y3=y
      sol3=sol
      sol3(i)=sol3(i)+serr(i,2)
      call lcmodel(nbodies,npt,tol,y3,sol3,time,itime,ans)
      m=m3 !store mass values to module
      chi2=0
      do j=1,n
         chi2=chi2+(flux(j)-ans(j))*(flux(j)-ans(j))/(ferr(j)*ferr(j))
      enddo
      dchi=abs(chi2-chi1)
   enddo

enddo
do i=1,nbodies
   isol=8+i
   soldiff(9)=max(1.0e-6,0.01*sol(isol))
   if(serr(isol,2).ne.0.0d0)then
      h=soldiff(9)
      serr(isol,2)=1.0d0/abs(dfridr(func,sol(isol),h,drerr))
   endif
   write(0,500) sol(isol),serr(isol,2)
enddo
isol=0
do i=1,nbodies
!   write(0,*) "M:",i,m(i),mserr(i,2)
   im=i
   mdiff=max(1.0e-7,0.01*m(im))
   if(mserr(i,2).ne.0.0d0)then
      h=mdiff
      mserr(i,2)=1.0d0/abs(dfridr(func,m(im),h,drerr))
   endif
   write(0,500) m(im),mserr(im,2)
enddo
im=0
allocate(ydiff(6))

do i=1,nbodies
   rstar=(m(1)/(4.0d0/3.0d0*Pi*sol(1)*1000.0d0)*MU)**(1.0d0/3.0d0)
   if(i.gt.1)then
      length=  (y(6*i-5)-y(1))*(y(6*i-5)-y(1))+ &
               (y(6*i-4)-y(2))*(y(6*i-4)-y(2))+ &
               (y(6*i-3)-y(3))*(y(6*i-3)-y(3))
      length=rstar/sqrt(LU2*length)/10.0d0 !1/10 the scaled Sun
!      write(0,*) "scale:",rstar,length
      ydiff(1)=length
      ydiff(2)=length
      ydiff(3)=length
      ydiff(4)=0.01*TU/LU !cm/s
      ydiff(5)=0.01*TU/LU
      ydiff(6)=0.01*TU/LU
   else
      ydiff(1)=rstar/1000.0d0
      ydiff(2)=rstar/1000.0d0
      ydiff(3)=rstar/1000.0d0
      ydiff(4)=0.01*TU/LU
      ydiff(5)=0.01*TU/LU
      ydiff(6)=0.01*TU/LU
   endif

   do j=1,6
      iy=6*(i-1)+j
      if(yserr(iy,2).ne.0.0d0)then
!         write(0,*) iy,y(iy),ydiff(j)
         h=ydiff(j)
         yserr(iy,2)=1.0d0/abs(dfridr(func,y(iy),h,drerr))
      endif
      write(0,500) y(iy),yserr(iy,2),ydiff(j)
   enddo
enddo
iy=0

write(0,*) "Exporting fit"
nplanet=nbodies-1
call exportfit(nplanet,sol,serr,err,y,yerr,yserr)

end program mcmcsetup7

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function func(x)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittingmod
implicit none
integer :: i,nunit
real(double) :: func,x
real(double), allocatable, dimension(:) :: ans,y3,sol3,m3

interface
   subroutine lcmodel(nbodies,npt,tol,y,sol,time,itime,ans)
      use precision
      implicit none
      integer, intent(in) :: nbodies,npt
      real(double), intent(in) :: tol
      real(double), dimension(:), intent(inout) :: y,ans
      real(double), dimension(:), intent(in) :: sol,time,itime
   end subroutine lcmodel
end interface

allocate(y3(nbodies2*6),sol3(8+nbodies2),m3(nbodies2))
do i=1,nbodies2*6
   y3(i)=y2(i)
enddo
do i=1,nbodies2
   m3(i)=m(i)
enddo
do i=1,8+nbodies2
   sol3(i)=sol2(i)
enddo

if(isol.gt.0) sol3(isol)=x
if(im.gt.0) m(im)=x
if(iy.gt.0) y3(iy)=x

!write(0,500) (m(i),i=1,nbodies2)
500  format(28(1PE10.3,1X))

allocate(ans(npt2))
call lcmodel(nbodies2,npt2,tol2,y3,sol3,time2,itime2,ans)

!nunit=10
!open(unit=nunit,file="junk.dat")
!do i=1,npt2
!   write(nunit,501) time2(i),flux2(i),ans(i)
!enddo
!close(nunit)
!501 format(5(1X,1PE17.10))
!write(0,*) "lcmodel in FCN done"
!read(5,*)

func=0.0d0
do i=1,npt2
   func=func+(ans(i)-flux2(i))*(ans(i)-flux2(i))/(ferr2(i)*ferr2(i))
enddo

do i=1,nbodies2
   m(i)=m3(i)
enddo

!write(0,*) x,func

return
end function func
