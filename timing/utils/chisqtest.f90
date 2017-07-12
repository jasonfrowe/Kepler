subroutine chisqtest(ntt,rtt,rtterr,chisq)
use precision
implicit none
integer :: ntt,i
real(double) :: chisq
real(double), dimension(:) :: rtt,rtterr

chisq=0.0d0
do i=1,ntt
   chisq=chisq+rtt(i)*rtt(i)/(rtterr(i)*rtterr(i))
!   write(0,*) chisq,rtt(i),rtterr(i)
!   read(5,*)
enddo

return
end subroutine chisqtest
