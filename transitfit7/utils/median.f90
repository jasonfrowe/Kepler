function median(n,data)
!Requires RQSORT (double precision) 
use precision
implicit none
!Import Vars
integer :: n
real(double) :: median
real(double), dimension(n) :: data
!Local Vars
integer, allocatable, dimension(:) :: p

allocate(p(n))

call rqsort(n,data,p)

if(n.le.0)then
	write(0,*) "No data for median"
	median=0.0
else
	median=data(p((n+1)/2))
endif

return
end function median