subroutine exportsol(nstar,Id,Ic,xcoo,ycoo,ngsol,gsol)
use precision
implicit none
integer :: nstar,ngsol,Id(nstar)
real(double) :: Ic(nstar),xcoo(nstar),ycoo(nstar),gsol(ngsol)
!local vars
integer :: i

write(6,500) (gsol(i),i=1,ngsol)
500 format(4(1PE17.10,1X))
do i=1,nstar
   write(6,501) Id(i),Ic(i),xcoo(i),ycoo(i)
enddo
501 format(I9,1X,3(1PE17.10,1X))

return
end subroutine exportsol
