subroutine fittransitmodel(nbodies,npt,tol,y,yserr,sol,serr,time,flux,ferr,itime)
use precision
use fittingmod
implicit none
integer npt,n,i,info,lwa,j
integer, target :: nbodies
integer, allocatable, dimension(:) :: iwa
real(double), target :: tol,tollm
real(double), dimension(:), target :: y,sol,time,flux,ferr,itime
real(double), dimension(:,:), target :: serr,yserr
real(double), allocatable, dimension(:) :: solfit,wa,fvec
external fcn

allocate(solfit(8+nbodies+nbodies+nbodies*6))

!how many variables are we fitting?
n=0
do i=1,8+nbodies
   if(serr(i,2).ne.0.0d0)then
      n=n+1
      solfit(n)=sol(i)
   endif
enddo
do i=1,nbodies
   if(mserr(i,2).ne.0.0d0)then
      n=n+1
      solfit(n)=m(i)
   endif
enddo
do i=1,nbodies*6
   if(yserr(i,2).ne.0.0d0)then
      n=n+1
      solfit(n)=y(i)
   endif
enddo

write(0,*) "n:",n

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

lwa=npt*n+5*npt*n
allocate(fvec(npt),iwa(n),wa(lwa))
tollm=1.0d-12
!!call lmdif1(fcn,m,n,x,fvec,tollm,info,iwa,wa,lwa)
write(0,*) "Calling lmdif1"
call lmdif1(fcn,npt,n,solfit,fvec,tollm,info,iwa,wa,lwa)
write(0,*) "info: ",info

n=0
do i=1,8+nbodies
   if(serr(i,2).ne.0.0d0)then
      n=n+1
      sol(i)=solfit(n)
      write(0,*) "out:",sol(i),solfit(n)
   endif
enddo
do i=1,nbodies
   if(mserr(i,2).ne.0.0d0)then
      n=n+1
      m(i)=solfit(n)
      write(0,*) "out:",m(i),solfit(n)
   endif
enddo
do i=1,nbodies*6
   if(yserr(i,2).ne.0.0d0)then
      n=n+1
      y(i)=solfit(n)
      write(0,*) "out:",y(i),solfit(n)
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
   subroutine lcmodel(nbodies,npt,tol,y,sol,time,itime,ans)
      use precision
      implicit none
      integer, intent(in) :: nbodies,npt
      real(double), intent(in) :: tol
      real(double), dimension(:), intent(inout) :: y,ans
      real(double), dimension(:), intent(in) :: sol,time,itime
   end subroutine lcmodel
end interface

!write(0,*) "FCN1"

allocate(sol3(8+nbodies2),y3(nbodies2*6))
allocate(m3(nbodies2))

j=0
do i=1,8+nbodies2
   sol3(i)=sol2(i)
enddo
do i=1,8+nbodies2
   if(serr2(i,2).ne.0.0d0)then
      j=j+1
      sol3(i)=x(j)
   endif
   write(0,501) sol3(i),sol2(i),sol3(i)-sol2(i),serr2(i,2)
enddo
do i=1,nbodies2 !save a copy of m? (mass is a global variable)
   m3(i)=m(i)
enddo
do i=1,nbodies2
   if(mserr(i,2).ne.0.0d0)then
      j=j+1
      m(i)=x(j)
   endif
   write(0,501) m(i),m3(i),m(i)-m3(i),mserr(i,2)
enddo
do i=1,nbodies2*6
   y3(i)=y2(i)
enddo
do i=1,nbodies2*6
   if(yserr2(i,2).ne.0.0d0)then
      j=j+1
      y3(i)=x(j)
   endif
   write(0,501) y3(i),y2(i),y3(i)-y2(i),yserr2(i,2)
enddo

!write(0,*) "lcmodel",nbodies,npt
call lcmodel(nbodies2,npt,tol2,y3,sol3,time2,itime2,fvec)
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
   fvec(i)=(fvec(i)-flux2(i))/ferr2(i)*yp
enddo

do i=1,nbodies2
   m(i)=m3(i) !restore m  (is this even necessary?)
enddo

return
end


