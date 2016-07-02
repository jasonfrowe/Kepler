subroutine findnei(nstar,xcoo,ycoo,neimax,numnei,nei,radnei)
use precision
implicit none
integer :: nstar,neimax
integer, dimension(:) :: numnei
integer, dimension(:,:) :: nei
real(double) :: radnei
real(double), dimension(:) :: xcoo,ycoo
!Local vars
integer i,j
real(double) :: dx,dy,radneid2

!precompute radnei/2
radneid2=radnei/2.0d0

!initialize neighbour count
numnei=0

do i=1,nstar
   do j=1,nstar
      if(i.ne.j)then
         dx=abs(xcoo(i)-xcoo(j))
         dy=abs(ycoo(i)-ycoo(j))
         if((dx.lt.radneid2).and.(dy.lt.radneid2))then
            numnei(i)=numnei(i)+1
            if(numnei(i).gt.neimax)then
               write(0,*) "neimax is too small!!"
               stop
            endif
            nei(i,numnei(i))=j
!            write(0,*) i,j,dx,dy
!            read(5,*)
         endif
      endif
   enddo
enddo

return
end subroutine findnei
