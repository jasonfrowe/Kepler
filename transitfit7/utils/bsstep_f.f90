!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine f(t,y,ydot)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: neq,npt,i,j
integer, parameter :: nmax=1000
real(double), allocatable, dimension(:) :: m
real(double) :: t,y(nmax*6),ydot(nmax*6),tmp(3),R

npt=nbodies2 !number of bodies from pointer 
allocate(m(npt))
m(1:npt)=m2(1:npt) !masses from pointer

do i=1,npt
!Okay.  for x,y,z the derivatives are simply vx,vy,vz.  simple.
   ydot(6*i-5)=y(6*i-2)
   ydot(6*i-4)=y(6*i-1)
   ydot(6*i-3)=y(6*i)
!Initialize tmp array
   do j=1,3
      tmp(j)=0.0d0
   enddo
   do j=1,npt
      if((j.ne.i).and.(m(j).gt.0.0d0))then
         R=(y(6*i-5)-y(6*j-5))*(y(6*i-5)-y(6*j-5))+ &
            (y(6*i-4)-y(6*j-4))*(y(6*i-4)-y(6*j-4))+ &
            (y(6*i-3)-y(6*j-3))*(y(6*i-3)-y(6*j-3))
         R=R**(3.0d0/2.0d0)
         tmp(1)=tmp(1)+G*m(j)*(y(6*j-5)-y(6*i-5))/R
         tmp(2)=tmp(2)+G*m(j)*(y(6*j-4)-y(6*i-4))/R
         tmp(3)=tmp(3)+G*m(j)*(y(6*j-3)-y(6*i-3))/R
      endif
   enddo
   ydot(6*i-2)=tmp(1)
   ydot(6*i-1)=tmp(2)
   ydot(6*i)=tmp(3)
enddo
 
return
end
