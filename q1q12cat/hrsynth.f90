program hrsynth
use precision
implicit none
integer :: iargc
character(80) :: catfile

!need at least 1 commandline argument
if(iargc().lt.1) then
   write(0,*) "Usage: hrsynth Q1Q12catalogue"
   stop
endif

!read in name of catalogue file from commandline
call getarg(1,catfile)



end program hrsynth
