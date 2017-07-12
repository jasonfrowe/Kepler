subroutine calcnbodies(inputsol,nbodies)
use precision
implicit none
integer :: nbodies,nunit,filestatus,rnbody
character(80) :: inputsol,command

nunit=10
nbodies=0
open(unit=nunit,file=inputsol,iostat=filestatus,status='old')

if(filestatus>0)then
   write(0,*) "Cannot open ",inputsol
   stop
endif

do
   read(nunit,*,iostat=filestatus) command
   if(filestatus == 0) then
      read(command(3:3),*,iostat=filestatus) rnbody
      if((filestatus == 0).and.(command(1:2).ne.'NL'))then
!         write(0,*) nbodies,rnbody
         if(rnbody.gt.nbodies) nbodies=rnbody
      endif
      cycle
   elseif(filestatus == -1) then
      exit
   else
      write(0,*) "File Error!!"
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo

close(nunit)

return
end subroutine
