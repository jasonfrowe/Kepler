program pickcandidates
use precision
implicit none
integer :: iargc,nunit,filestatus,ncol,i,npt,cflag
real(double), allocatable, dimension(:) :: values
character(80) :: pickfile

interface
   subroutine checks(ncol,values,cflag)
      use precision
      implicit none
      integer, intent(inout) :: ncol,cflag
      real(double), dimension(:), intent(inout) :: values
   end subroutine checks
end interface

if(iargc().lt.1) then  !check that we have sufficient commandline arguments
   write(0,*) "Usage: pickcandidates lptarg_s1_pt16_sn7.dat>"
   stop
endif

call getarg(1,pickfile)  !get filename for photometry

nunit=10 !unit number for data spectrum
open(unit=nunit,file=pickfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",pickfile
   stop
endif

ncol=24
allocate(values(ncol))

npt=0 !count number of candidates
!read through the list of candidates and pick the ones that meet our criteria
do
   read(nunit,*,iostat=filestatus) (values(i),i=1,ncol)
   if(filestatus == 0) then
      npt=npt+1
      call checks(ncol,values,cflag)
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

end program pickcandidate

!CCCCCCCCCCCCCCCC
subroutine checks(ncol,values,cflag)
!CCCCCCCCCCCCCCCC
use precision
implicit none
integer :: ncol, cflag
real(double), dimension(:) :: values

if(values(8) > 7.0)then !check S/N
   if((values(10).gt.0.4995).or.(values(10).gt.0.5))then

   else
      cflag=2 !failed fbmax1 test

else
   cflag=1 !failed S/N
endif

return
end
