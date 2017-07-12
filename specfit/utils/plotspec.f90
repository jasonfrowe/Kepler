!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine plotspec(npt,wv,flux)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer npt,i
real(double), dimension(:) :: wv,flux
real, allocatable, dimension(:) :: xp,yp,rj

allocate(rj(4),xp(npt),yp(npt))

do i=1,npt
   xp(i)=real(wv(i))
   yp(i)=real(flux(i))
enddo
rj(1)=minval(xp)
rj(2)=maxval(xp)
rj(3)=minval(yp)-0.1*(maxval(yp)-minval(yp))
rj(4)=maxval(yp)+0.1*(maxval(yp)-minval(yp))

call pgwindow(rj(1),rj(2),rj(3),rj(4))
!write(0,*) rj
call pgbox("BCNTS1",0.0,0,"BCNTS1",0.0,0)
call pglabel("Wavelength (\A)","Flux","")

call pgline(npt,xp,yp)

deallocate(rj,xp,yp)

return
end subroutine plotspec
