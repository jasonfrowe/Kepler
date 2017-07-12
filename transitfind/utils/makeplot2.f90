!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine makeplot(nplot,nstep,freqs,p,obsfile,epo,bper,bpow,depth,    &
   pmean,std,qtran,npt,time,flux)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nplot,nstep,npt
real(double) :: epo,bper,phase,bpow,depth,pmean,std,qtran
real(double), dimension(:) :: freqs,p,time,flux
real(double), allocatable, dimension(:) :: pers
character(80) :: obsfile,txtout

select case(nplot)
   case(1)
      call pgopen('/xserve')
   case(2)
      call pgopen('find.ps/vcps')
   case(3)
      call pgopen('find.png/png')
end select

call pgsch(2.9)
call pgsubp(1,4)
call pgvport(0.1,0.95,0.2,0.8) !gives enough room for labels

call pgpage()

if(npt.ge.3)then

   allocate(pers(nstep))
   pers=1.0d0/freqs
   call plot(nstep,pers,p,obsfile,"Period (d)","P",0,0.0d0,0.0d0)
   call pgpage()
   call plot(npt,time,flux,obsfile,"BJD-2454900","Flux",1,epo,bper)
   call pgpage()
   phase=epo/bper-int(epo/bper)
   if(phase.lt.0.0d0) phase=phase+1.0d0
   call plotph(npt,time,flux,bper,phase)
   call pgpage()
   call pgwindow(0.0,1.0,0.0,1.0)

   write(txtout,501) "Per: ",bper
   501 format(A5,1X,1PE17.10)
   call pgptxt(0.1,0.9,0.0,0.5,txtout)
   write(txtout,501) "Epo: ",epo
   call pgptxt(0.1,0.77,0.0,0.5,txtout)
   write(txtout,501) "Pow: ",bpow
   call pgptxt(0.1,0.64,0.0,0.5,txtout)
   write(txtout,501) "S/N: ",(depth-pmean)/std*sqrt(qtran*dble(npt))
   call pgptxt(0.1,0.51,0.0,0.5,txtout)
   write(txtout,501) "Dur: ",qtran*bper
   call pgptxt(0.1,0.38,0.0,0.5,txtout)
   write(txtout,501) "Dep: ",depth*1.0e6
   call pgptxt(0.1,0.25,0.0,0.5,txtout)

   call plottrans(npt,time,flux,bper,phase,qtran)

endif

call pgclos()

return
end
