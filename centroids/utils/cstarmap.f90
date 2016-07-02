subroutine cstarmap(nstar,xcoo,ycoo,naxes,Image,starmap,bpix)
use precision
implicit none
integer :: nstar
integer, dimension(2) :: naxes
integer, dimension(:,:) :: starmap
real(double) :: bpix
real(double), dimension(:) :: xcoo, ycoo
real(double), dimension(:,:) :: Image
!local vars
integer :: i,j,k
real(double) :: d,di,dj,dmin,dmininit

!a big value for the distance
dmininit=(dble(max(naxes(1),naxes(2))))**2.0

do i=1,naxes(1)
   di=dble(i)
   do j=1,naxes(2)
      dj=dble(j)
      if(Image(i,j).lt.bpix)then
         dmin=dmininit
         do k=1,nstar
            d= sqrt((xcoo(k)-di)**2.0 + (ycoo(k)-dj)**2.0)
            if(d.lt.dmin)then
               dmin=d
               starmap(i,j)=k
            endif
         enddo
!         write(0,*) i,j,starmap(i,j),image(i,j)
!         read(5,*)
      endif
   enddo
enddo

return
end subroutine cstarmap
