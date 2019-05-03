module ocmod
   use precision, only: double
   implicit none
   integer :: ntmidmax
   integer, allocatable, dimension(:) :: ntmid
   real(double), allocatable, dimension(:,:) :: tmid
end module ocmod
