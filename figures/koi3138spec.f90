program koi3138spec
use precision
implicit none
integer :: nmax,iargc,nunit,filestatus,nptR,i,nptB
real, allocatable, dimension(:) :: lamR,fluxR,lamB,fluxB,rj,lx,ly,px,py
character(80) :: redspec, bluespec

if(iargc().lt.2)then
   write(0,*) "Usage: koi3138srad redspectrum bluespectrum"
   write(0,*) " redspectrum - txt file with red spectrum"
   write(0,*) " bluespectrum - txt file with blue spectrum"
   stop
endif

nmax=5000

call getarg(1,redspec)
allocate(lamR(nmax),fluxR(nmax))

nunit=10 !unit number for data spectrum
open(unit=nunit,file=redspec,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",redspec
   stop
endif
i=1
do
   if(i.gt.nmax)then !check that we have enough allocated memory
      write(0,*) "Critical Error: Increase nmax"
      stop
   endif
   read(nunit,*,iostat=filestatus) lamR(i),fluxR(i)
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
close(nunit)
nptR=i-1  !number of chains read from stellar MCMC work
write(0,*) "nptR: ",nptR
fluxR=fluxR*1e16 !scale up fluxes

call getarg(2,bluespec)
allocate(lamB(nmax),fluxB(nmax))
open(unit=nunit,file=bluespec,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",bluespec
   stop
endif
i=1
do
   if(i.gt.nmax)then !check that we have enough allocated memory
      write(0,*) "Critical Error: Increase nmax"
      stop
   endif
   read(nunit,*,iostat=filestatus) lamB(i),fluxB(i)
   if(filestatus == 0) then
      i=i+1
      cycle
   elseif(filestatus == -1) then
      exit
   else
      write(0,*) "File Error!!"
      write(0,900) "iostat: ",filestatus
      stop
   endif
enddo
close(nunit)
nptB=i-1  !number of chains read from stellar MCMC work
write(0,*) "nptB: ",nptB
close(nunit)
fluxB=fluxB*1e16 !scale up fluxes

call pgopen('?') !open plotting device
call pgpage()
call pgask(.false.)
call PGPAP ( 8.0 ,1.0) !use a square 8" across
call pgsubp(1,2)
call pgpage()
call pgsch(1.5) !make the font a bit bigger
call pgslw(3)  !make the lines a bit thicker
call pgvport(0.15,0.85,0.15,0.95) !make room around the edges for labels

allocate(rj(4))
rj(1)=minval(lamR(1:nptR)) !scale for plot
rj(2)=maxval(lamR(1:nptR))
rj(3)=minval(fluxR(1:nptR))-0.1*(maxval(fluxR(1:nptR))-minval(fluxR(1:nptR)))
rj(4)=maxval(fluxR(1:nptR))+0.1*(maxval(fluxR(1:nptR))-minval(fluxR(1:nptR)))
call pgwindow(rj(1),rj(2),rj(3),rj(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS1",0.0,0)
call pglabel("Wavelength (\A)","Flux (10\u-16\d erg s\u-1\d cm\u-2\d \A\u-1\d)","")
allocate(px(4),py(4))
px(1)=9200.0
py(1)=5.0!rj(3)
px(2)=rj(2)
py(2)=5.0!rj(3)
px(3)=rj(2)
py(3)=rj(4)
px(4)=9200.0
py(4)=rj(4)
call pgsci(2)
call PGSFS(3)
call PGSHS(45.0, 2.0, 0.0)
call pgpoly(4,px,py)
call PGSHS(45.0, 1.0, 0.0)
call pgsci(1)
call PGSFS(1)

call pgline(nptR,lamR,fluxR)

allocate(lx(2),ly(2))
call pgsch(1.0)
call pgslw(2)
call pgptxt(7682.0,1.0,0.0,0.5,"\frKI\fn")
lx(1)=7665.0
lx(2)=7665.0
ly(1)=1.5
ly(2)=1.8
call pgline(2,lx,ly)
lx(1)=7699.0
lx(2)=7699.0
call pgline(2,lx,ly)
lx(1)=7665.0
lx(2)=7699.0
ly(1)=1.5
ly(2)=1.5
call pgline(2,lx,ly)
call pgptxt(8189.0,3.0,0.0,0.5,"\frNaI\fn")
lx(1)=8183.0
lx(2)=8183.0
ly(1)=3.5
ly(2)=3.8
call pgline(2,lx,ly)
lx(1)=8195.0
lx(2)=8195.0
call pgline(2,lx,ly)
lx(1)=8183.0
lx(2)=8195.0
ly(1)=3.5
ly(2)=3.5
call pgline(2,lx,ly)
call pgptxt(7128.0,6.8,0.0,0.5,"\frTiO\fn")
lx(1)=7053.0
lx(2)=7203.0
ly(1)=6.5
ly(2)=6.5
call pgline(2,lx,ly)
call pgptxt(7741.0,9.0,0.0,0.5,"\frTiO\fn")
lx(1)=7666.0
lx(2)=7816.0
ly(1)=8.7
ly(2)=8.7
call pgline(2,lx,ly)
call pgptxt(8356.0,9.5,0.0,0.5,"\frTiO\fn")
lx(1)=8206.0
lx(2)=8506.0
ly(1)=9.2
ly(2)=9.2
call pgline(2,lx,ly)
call pgptxt(7417.0,8.3,0.0,0.5,"\frVO\fn")
lx(1)=7334.0
lx(2)=7500.0
ly(1)=8.0
ly(2)=8.0
call pgline(2,lx,ly)
call pgptxt(7920.0,8.3,0.0,0.5,"\frVO\fn")
lx(1)=7851.0
lx(2)=7990.0
ly(1)=8.0
ly(2)=8.0
call pgline(2,lx,ly)
call pgptxt(8580.0,4.0,0.0,0.5,"\frCaII\fn")
lx(1)=8498.03
lx(2)=8662.14
ly(1)=4.5
ly(2)=4.5
call pgline(2,lx,ly)
lx(1)=8498.03
lx(2)=8498.03
ly(1)=4.5
ly(2)=4.8
call pgline(2,lx,ly)
lx(1)=8542.09
lx(2)=8542.09
ly(1)=4.5
ly(2)=4.8
call pgline(2,lx,ly)
lx(1)=8662.14
lx(2)=8662.14
ly(1)=4.5
ly(2)=4.8
call pgline(2,lx,ly)

call pgpage()
call pgsch(1.5)
call pgslw(3)
rj(1)=minval(lamB(1:nptB)) !scale for plot
rj(2)=maxval(lamB(1:nptB))
rj(3)=minval(fluxB(1:nptB))-0.1*(maxval(fluxB(1:nptB))-minval(fluxB(1:nptB)))
rj(4)=maxval(fluxB(1:nptB))+0.1*(maxval(fluxB(1:nptB))-minval(fluxB(1:nptB)))
call pgwindow(rj(1),rj(2),rj(3),rj(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS1",0.0,0)
call pglabel("Wavelength (\A)","Flux (10\u-16\d erg s\u-1\d cm\u-2\d \A\u-1\d)","")
call pgline(nptB,lamB,fluxB)

call pgsch(1.0)
call pgslw(2)
call pgptxt(4227.0,1.5,0.0,0.5,"\frCaI\fn")
lx(1)=4227.0
lx(2)=4227.0
ly(1)=1.2
ly(2)=0.9
call pgline(2,lx,ly)
call pgptxt(4773.0,2.5,0.0,0.5,"\frTiO\fn")
lx(1)=4750.0
lx(2)=4790.0
ly(1)=2.2
ly(2)=2.2
call pgline(2,lx,ly)
call pgptxt(5010.0,3.5,0.0,0.5,"\frTiO\fn")
lx(1)=4950.0
lx(2)=5070.0
ly(1)=3.2
ly(2)=3.2
call pgline(2,lx,ly)
call pgptxt(5185.0,3.5,0.0,0.5,"\frTiO\fn")
lx(1)=5160.0
lx(2)=5210.0
ly(1)=3.2
ly(2)=3.2
call pgline(2,lx,ly)
call pgptxt(5467.0,4.3,0.0,0.5,"\frTiO\fn")
lx(1)=5445.0
lx(2)=5490.0
ly(1)=4.0
ly(2)=4.0
call pgline(2,lx,ly)
call pgptxt(5530.0,5.2,0.0,0.5,"\frCaOH\fn")
lx(1)=5500.0
lx(2)=5560.0
ly(1)=4.9
ly(2)=4.9
call pgline(2,lx,ly)
call pgptxt(5892.9375,0.0,0.0,0.5,"\frNaD\fn")
lx(1)=5889.951
lx(2)=5889.951
ly(1)=0.3
ly(2)=0.6
call pgline(2,lx,ly)
lx(1)=5895.924
lx(2)=5895.924
ly(1)=0.3
ly(2)=0.6
call pgline(2,lx,ly)
lx(1)=5889.951
lx(2)=5895.924
ly(1)=0.3
ly(2)=0.3
call pgline(2,lx,ly)

call pgclos()

end program koi3138spec
