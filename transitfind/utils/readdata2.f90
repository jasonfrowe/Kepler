!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nmax,npt
real(double) :: ztime
real(double), dimension(:) :: time,flux,ferr,itime
character(80) :: obsfile
!local vars
integer :: nunit,filestatus,i
real(double) :: t,f,e,sec2day

sec2day=86400.0d0

nunit=10
open(unit=nunit,file=obsfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",obsfile
   stop
endif

i=0
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) t,f,e
   if(filestatus == 0) then
      i=i+1
      time(i)=t-ztime+0.5d0
      flux(i)=f+1.0
      ferr(i)=e
      itime(i)=1765.5/sec2day
      !itime(i)=58.85/sec2day !short cadence
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i

return
end subroutine readdata2
