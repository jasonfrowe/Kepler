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
real(double) :: t,f,e,sec2day,it,min2day
character(400) :: line
character(20) :: zeros

sec2day=86400.0d0 !number of seconds in a day
min2day=24.0*60.0 !number of minutes in a day
write(zeros, '(*(I2))') (0, i=1, 4) !write a bunch of zeros  

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
   read(nunit,'(A)',iostat=filestatus) line
   !write(0,*) t,f,e,it
   !it=0.0 !Forces the use of Long Cadence
   if(filestatus == 0) then
      line=trim(line) // zeros  !trick to always read in at least 4 columns. 
      read(line,*,iostat=filestatus) t,f,e,it
      if(filestatus == 0) then
         i=i+1
         time(i)=t-ztime+0.5d0
         flux(i)=f+1.0
         ferr(i)=e
         if ((it.lt.1.0d-7).and.(it.gt.-1.0d-7)) then
            itime(i)=1765.5/sec2day !long cadence
         elseif (it.lt.0.0) then
            itime(i)=58.85/sec2day !short cadence
         else
            itime(i)=it !Units are now in days   !/min2day !itime is expected to be in minutes when read in, convert to days.
         endif
         !write(0,500) i,time(i),flux(i),ferr(i),itime(i)
         500 format(I6,1X,F10.5,1X,F8.6,1X,F8.6,1X,F8.6)
         else
            write(0,*) "File Error!! Line:",i+1
            write(0,900) "iostat: ",filestatus
            stop
         endif
      !read(5,*)
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
