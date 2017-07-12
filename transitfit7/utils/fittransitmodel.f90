subroutine fittransitmodel(nbodies,npt,tol,sol,serr,time,flux,ferr,itime)
use precision
use fittingmod
implicit none
integer npt,n,i,info,lwa,j
integer, target :: nbodies
integer, allocatable, dimension(:) :: iwa
real(double), target :: tol,tollm
real(double), dimension(:), target :: sol,time,flux,ferr,itime
real(double), dimension(:,:), target :: serr
real(double), allocatable, dimension(:) :: solfit,wa,fvec
external fcn

allocate(solfit(7+nbodies*7))

!how many variables are we fitting?
n=0
do i=1,7+nbodies*7
   if(serr(i,2).ne.0.0d0)then
      n=n+1
      solfit(n)=sol(i)
   endif
enddo
write(0,*) "n:",n

!assign pointers from module that is shared with FCN
tol2 => tol
nbodies2 => nbodies
sol2 => sol
serr2 => serr
time2 => time
flux2 => flux
ferr2 => ferr
itime2 => itime

lwa=npt*n+5*npt*n
allocate(fvec(npt),iwa(n),wa(lwa))
tollm=1.0d-12
!!call lmdif1(fcn,m,n,x,fvec,tollm,info,iwa,wa,lwa)
write(0,*) "Calling lmdif1"
call lmdif1(fcn,npt,n,solfit,fvec,tollm,info,iwa,wa,lwa)
write(0,*) "info: ",info

n=0
do i=1,7+nbodies*7
   if(serr(i,2).ne.0.0d0)then
      n=n+1
      sol(i)=solfit(n)
      write(0,*) "out:",sol(i),solfit(n)
   endif
enddo

return
end subroutine fittransitmodel

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn(npt,n,x,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittingmod
implicit none
integer :: npt,n,iflag,i,j,nunit
real(double) :: yp
real(double), allocatable, dimension(:) :: sol3,m3,y3
real(double), dimension(n) :: x
real(double), dimension(npt) :: fvec

interface
   subroutine lcmodel(nbodies,npt,tol,sol,time,itime,ans)
      use precision
      implicit none
      integer, intent(in) :: nbodies,npt
      real(double), intent(in) :: tol
      real(double), dimension(:), intent(inout) :: ans
      real(double), dimension(:), intent(in) :: sol,time,itime
   end subroutine lcmodel
end interface

!write(0,*) "FCN1"

allocate(sol3(7+nbodies2*7))

j=0
do i=1,7+nbodies2*7
   sol3(i)=sol2(i)
enddo
do i=1,7+nbodies2*7
   if(serr2(i,2).ne.0.0d0)then
      j=j+1
      sol3(i)=x(j)
   endif
   write(0,501) sol3(i),sol2(i),sol3(i)-sol2(i),serr2(i,2)
enddo

!write(0,*) "lcmodel",nbodies,npt
call lcmodel(nbodies2,npt,tol2,sol3,time2,itime2,fvec)
!write(0,*) "lcmodel done"

nunit=10
open(unit=nunit,file="junk.dat")
do i=1,npt
   write(nunit,501) time2(i),flux2(i),fvec(i)
enddo
close(nunit)
501 format(5(1X,1PE17.10))
write(0,*) "lcmodel in FCN done"
!read(5,*)

yp=1.0d0
do i=1,npt
   fvec(i)=(fvec(i)-flux2(i))/ferr2(i)*yp  !nope.. this is wrong..
enddo

return
end


