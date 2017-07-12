program k10phaseplot
use precision
implicit none
integer :: nmax,npt,iargc,nplot,i,nbin,ibin
integer, allocatable, dimension(:) :: pn
real :: tr
real, allocatable, dimension(:) :: bb,px,py,pe
real(double) :: P,T0,phase,dnbin,err2,zpt
real(double), allocatable, dimension(:) :: time,flux,ferr,tsum,fsum,    &
   esum,tmodel,t1,t2
character(80) :: filename

P=8.3749120981d-01
T0=6.4574763507d+01
zpt=-3.2389260896E-06

if(iargc().lt.1)then
   write(0,*) "Usage: k10phaseplot <lightcurve>"
   write(0,*) " <lightcurve> : 3 column file with time,flux,err"
   stop
endif

!get commandland argument
call getarg(1,filename)

nmax=2000000
allocate(time(nmax),flux(nmax),ferr(nmax),tmodel(nmax))
call readphot(filename,nmax,npt,time,flux,ferr,tmodel)
write(0,*) "Number of points read: ",npt

time=time-t0
flux=flux+1.0
allocate(px(npt),py(npt))
do i=1,npt
   phase=time(i)/P-floor(time(i)/P)
   px(i)=real(phase)
   py(i)=real(flux(i))
enddo
nplot=npt

!get plot bounds
allocate(bb(4))
bb(1)=0.0d0!minval(px(1:nplot))
bb(2)=1.0d0!maxval(px(1:nplot))
bb(3)=minval(py(1:nplot))
bb(4)=maxval(py(1:nplot))
tr=bb(4)-bb(3)
bb(3)=bb(3)-0.1*tr
bb(4)=bb(4)+0.1*tr
write(0,*) bb

!open PGPLOT device
call pgopen('?')!('/xserve')  !'?' lets the user choose the device.
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgsubp(1,2)
call pgpage() !create a fresh page
call pgslw(3) !thicker lines
call pgsch(1.5) !bigger text

call pgvport(0.15,0.95,0.15,0.95) !make room around the edges for labels
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
call pglabel("Phase","Normalized Flux","")

call pgpt(npt,px,py,-1)

nbin=50
dnbin=dble(nbin)
allocate(tsum(nbin),fsum(nbin),esum(nbin))
!initialize arrays to zero
tsum=0.0d0
fsum=0.0d0
esum=0.0d0
do i=1,npt
   phase=time(i)/P-floor(time(i)/P)
   ibin=int(phase*dnbin)+1
   if((ibin.gt.0).and.(ibin.le.nbin))then
      err2=ferr(i)*ferr(i)
      tsum(ibin)=tsum(ibin)+phase/err2
      fsum(ibin)=fsum(ibin)+flux(i)/err2
      esum(ibin)=esum(ibin)+1.0d0/err2
   endif
enddo
do i=1,nbin
   if(esum(i).gt.0.0)then
      tsum(i)=tsum(i)/esum(i)
      fsum(i)=fsum(i)/esum(i)
      esum(i)=1.0d0/sqrt(esum(i))
!      write(0,*) tsum(i),fsum(i),esum(i)
   endif
enddo

call pgsci(2)
deallocate(px,py)
allocate(px(nbin),py(nbin),pe(nbin))
px(1:nbin)=real(tsum(1:nbin))
py(1:nbin)=real(fsum(1:nbin))
pe(1:nbin)=real(esum(1:nbin))
call pgpt(nbin,px,py,17)
call pgsci(1)

call pgpage()
bb(1)=0.0
bb(2)=1.0
bb(3)=-20.0
bb(4)= 20.0
call pgvport(0.15,0.95,0.15,0.95) !make room around the edges for labels
call pgwindow(bb(1),bb(2),bb(3),bb(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS",0.0,0)
call pglabel("Phase","Change in Flux","")
py=1.0e6*(py-1.0)

call pgpt(nbin,px,py,17)
call pgerrb(5,nbin,px,py,pe,1.0)

deallocate(px,py)
allocate(px(npt),py(npt),t1(npt),t2(npt))
do i=1,npt
   phase=time(i)/P-floor(time(i)/P)
   t1(i)=phase
   t2(i)=1.0d6*(tmodel(i)-1.0d0-zpt)
enddo
allocate(pn(npt))
call rqsort(npt,t1,pn)

do i=1,npt
   px(i)=real(t1(pn(i)))
   py(i)=real(t2(pn(i)))
enddo

call pgsci(2)
call pgline(npt,px,py)
call pgsci(1)

call pgclos()

end program k10phaseplot

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readphot(filename,nmax,npt,time,flux,ferr,tmodel)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: nmax,npt
real(double), dimension(nmax) :: time,flux,ferr,tmodel
character(80) :: filename
!local vars
integer :: nunit,filestatus,i

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

i=1
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) time(i),flux(i),ferr(i),tmodel(i)
   if(filestatus == 0) then
      i=i+1
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i-1

time=time-54900.0d0+0.5d0

return
end
