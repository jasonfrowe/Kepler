program kpixread
use precision
implicit none
integer :: iargc

if(iargc().lt.1)then
   write(0,*) "Usage: kpixread input.list"
   stop
endif

end program kpixread
