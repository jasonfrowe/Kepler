subroutine readstars(coofile,nstarmax,nstar,id,xcoo,ycoo)
use precision
implicit none
!import vars
integer :: nstarmax,nstar
integer, dimension(nstarmax) :: id
real(double), dimension(nstarmax) :: xcoo,ycoo
character(80) :: coofile
!local vars
integer :: nunit,filestatus,i,it
real(double) :: tx,ty

nunit=10
open(unit=nunit,file=coofile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",coofile
   stop
endif

i=0
do
   if(i.gt.nstarmax)then
      write(0,*) "Increase nstarmax to match data points"
      write(0,*) "nstarmax: ",nstarmax
      stop
   endif
   read(nunit,*,iostat=filestatus) it,tx,ty
   if(filestatus == 0) then
      i=i+1
      id(i)=it
      xcoo(i)=tx
      ycoo(i)=ty
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
nstar=i

return
end subroutine readstars
