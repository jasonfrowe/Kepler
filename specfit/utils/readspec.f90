subroutine readspec(nunit,nmax,npt,wv,flux,iflag)
use precision
implicit none
integer :: npt,iflag,nunit,i,nmax,filestatus
real(double), dimension(:) :: wv,flux

iflag=0
i=npt+1
do
   if(i.gt.nmax)then !ran out of array space
      npt=i-1
      iflag=1
      return
   endif
   read(nunit,*,iostat=filestatus) wv(i),flux(i)
   if(filestatus == 0) then
      i=i+1
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
npt=i-1 !number of colours read in

return
end subroutine readspec
