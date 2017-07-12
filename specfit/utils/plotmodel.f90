subroutine plotmodel(nmodel,wmod,fmod,minw,maxw)
use precision
implicit none
integer :: nmodel,i,nplot
real(double), dimension(:) :: wmod,fmod
real(double) :: minw,maxw
real, allocatable, dimension(:) :: rj,xp,yp

allocate(rj(4),xp(nmodel),yp(nmodel))

nplot=0
do i=1,nmodel
   if((wmod(i).gt.minw).and.(wmod(i).lt.maxw))then
      nplot=nplot+1
      xp(nplot)=real(wmod(i))
      yp(nplot)=real(fmod(i))
   endif
enddo

rj(1)=real(minw)
rj(2)=real(maxw)
rj(3)=minval(yp(1:nplot))-0.1*(maxval(yp(1:nplot))-minval(yp(1:nplot)))
rj(4)=maxval(yp(1:nplot))+0.1*(maxval(yp(1:nplot))-minval(yp(1:nplot)))
!write(0,*) "rj:",rj

call pgwindow(rj(1),rj(2),rj(3),rj(4))
!write(0,*) rj
call pgbox("BCNTS1",0.0,0,"BCNTS1",0.0,0)
call pglabel("Wavelength (\A)","Flux","")

call pgline(nplot,xp,yp)
deallocate(rj,xp,yp)

return
end subroutine plotmodel
