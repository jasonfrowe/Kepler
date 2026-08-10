function meddiff(npt,x)
use precision
implicit none
!input data
integer :: npt
real(double) :: meddiff
real(double), dimension(npt) :: x
!local data
integer :: i
integer, allocatable, dimension(:) :: p
real(double), allocatable, dimension(:) :: dd

allocate(p(npt-1),dd(npt-1))
do i=1,npt-1
        dd(i)=abs(x(i)-x(i+1))
enddo

call rqsort(npt-1,dd,p)
meddiff=dd(p((npt-1)/2))

return
end
