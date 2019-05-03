subroutine readkeplc(photfile,npt,time,flux,ferr,itime)
use precision
implicit none
integer :: npt,nunit,i,filestatus,nmax
real(double), dimension(:) :: time,flux,ferr,itime
real(double) :: Keplertime,sec2day
character(80) :: photfile

Keplertime=54900.0
sec2day=86400.0d0

nunit=10
open(unit=nunit,file=photfile,iostat=filestatus,status='old')

if(filestatus>0)then
   write(0,*) "Cannot open ",photfile
   stop
endif

i=1
do
!   write(0,*) i
   read(nunit,*,iostat=filestatus) time(i),flux(i),ferr(i)
   if(filestatus == 0) then
      time(i)=time(i)-Keplertime
      time(i)=time(i)+0.5d0 !MJD half day offset
      flux(i)=flux(i)+1.0!-2.5*log10(mag(i)+1.0d0)
      itime(i)=1765.5/sec2day
!      itime(i)=58.85/sec2day !short cadence
      i=i+1
   elseif(filestatus == -1) then
      exit
   else
      write(0,*) "File Error!! on line ",i
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
npt=i-1
close(nunit)

end subroutine
