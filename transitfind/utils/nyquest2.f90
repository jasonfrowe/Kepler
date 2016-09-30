!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine nyquest(npt,time,nyq)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: npt,ndt,i
integer, allocatable, dimension(:) :: p
real(double) :: nyq,mode
real(double), dimension(:) :: time
real(double), allocatable, dimension(:) :: dt

allocate(dt(npt),p(npt)) !allocate arrays

ndt=0 !initalize counter
do i=2,npt
   ndt=ndt+1
   dt(ndt)=time(i)-time(i-1) !difference between times
enddo

call rqsort(ndt,dt,p) !sort data
mode=dt(p(ndt/2)) !get median value

nyq=1.0/(2.0*mode) !calculate nyquest freq

!write(0,*) "nyq:",npt,nyq

return
end
