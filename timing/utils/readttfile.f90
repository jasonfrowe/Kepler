subroutine readttfile(ttfile,ntt,rtime,rtt,rtterr)
use precision
implicit none
integer :: ntt,nunit,filestatus,i,j
real(double) :: sigcut
real(double), dimension(:) :: rtime,rtt,rtterr
character(80) :: ttfile

nunit=11
open(unit=nunit,file=ttfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",ttfile
   stop
endif

ntt=1
do
   read(nunit,*,iostat=filestatus) rtime(ntt),rtt(ntt),rtterr(ntt)
   if(filestatus == 0) then
      if((rtterr(ntt).gt.0.0d0).and.(rtterr(ntt).lt.0.1)) ntt=ntt+1 !count number of files
      cycle
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!!",ttfile
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit)
ntt=ntt-1

sigcut=5.0
j=1
do i=2,ntt-1
   if((abs(rtt(i)-rtterr(i)).ge.sigcut).and.    &
    (abs(rtt(i-1)-rtterr(i-1)).lt.sigcut).and.  &
    (abs(rtt(i+1)-rtterr(i+1)).lt.sigcut))then
      continue
   else
      j=j+1
   endif
   rtime(j)=rtime(i)
   rtt(j)=rtt(i)
   rtterr(j)=rtterr(i)
enddo
j=j+1
rtime(j)=rtime(ntt)
rtt(j)=rtt(ntt)
rtterr(j)=rtterr(ntt)
ntt=j

close(nunit)

return
end subroutine readttfile
