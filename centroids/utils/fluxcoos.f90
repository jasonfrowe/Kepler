subroutine fluxcoos(nstar,xcoo,ycoo,naxes,Image,bpix)
use precision
implicit none
integer :: nstar
integer, dimension(2) :: naxes
real(double) :: bpix
real(double), dimension(:) :: xcoo,ycoo
real(double), dimension(:,:) :: Image
!local vars
integer :: i,j,k,dp,dj,dk
real(double) :: sumx,sumy,fsum,xc,yc

dp=2 !dp=2 produces a 5x5 window to find centroid center

do i=1,nstar !loop over all stars
   sumx=0.0d0 !initialize to zero
   sumy=0.0d0
   fsum=0.0d0
   do j=max(1,int(xcoo(i))-dp),min(naxes(1),int(xcoo(i))+dp)
      dj=dble(j)
      do k=max(1,int(ycoo(i))-dp),min(naxes(2),int(ycoo(i))+dp)
         dk=dble(k)
         if(Image(j,k).lt.bpix)then
            sumx=sumx+Image(j,k)*dj
            sumy=sumy+Image(j,k)*dk
            fsum=fsum+Image(j,k)
         endif
      enddo
   enddo
   xc=sumx/fsum
   yc=sumy/fsum

   xcoo(i)=xc
   ycoo(i)=yc
!   write(0,*) xcoo(i),xc
!   write(0,*) ycoo(i),yc
!   read(5,*)
enddo


return
end subroutine fluxcoos
