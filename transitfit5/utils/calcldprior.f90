function calcldprior(npriors,teff,logg,feh,ntype)
use precision
implicit none
!import vars
integer :: npriors,ntype
real(double) :: teff, logg, feh
real, dimension(:), pointer :: calcldprior
!local vars
integer :: i
real, dimension(:), allocatable, target :: ldprior

allocate(ldprior(npriors))

do i=1,npriors
	ldprior(i)=dble(i)
enddo

calcldprior => ldprior

return
end function calcldprior