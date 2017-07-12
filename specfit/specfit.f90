program specfit
use precision
implicit none
integer :: npt,nmax,iflag,iargc,nunit,filestatus,i,nmodelmax,nmodel,nfit
integer, allocatable, dimension(:) :: q
real(double) :: minw,maxw,norm,chisq,rv
real(double), allocatable, dimension(:) :: wv,flux,wv2,flux2,wmod,fmod, &
   wmod2,fmod2,fmodbin,ratio,diff
character(80) :: specfile,modelfile,cinput

!below are all the interfaces to allow dynamic arrays.
interface
   subroutine readspec(nunit,nmax,npt,wv,flux,iflag)
      use precision
      implicit none
      integer, intent(inout) :: npt,iflag,nunit,nmax
      real(double), dimension(:), intent(inout) :: wv,flux
   end subroutine readspec
end interface
interface
   subroutine plotspec(npt,wv,flux)
      use precision
      implicit none
      integer, intent(inout) :: npt
      real(double), dimension(:), intent(inout) :: wv, flux
   end subroutine plotspec
end interface
interface
   subroutine readmodel(nunit,nmodelmax,nmodel,wmod,fmod,iflag)
      use precision
      implicit none
      integer, intent(inout) :: nunit,nmodelmax,nmodel,iflag
      real(double), dimension(:), intent(inout) :: wmod, fmod
   end subroutine
end interface
interface
   subroutine plotmodel(nmodel,wmod,fmod,minw,maxw)
      use precision
      implicit none
      integer, intent(inout) :: nmodel
      real(double), dimension(:), intent(inout) :: wmod,fmod
      real(double), intent(inout) :: minw,maxw
   end subroutine
end interface
interface
   subroutine binmodel(npt,wv,nmodel,wmod,fmod,fmodbin,rv)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmodel
      real(double), intent(inout) :: rv
      real(double), dimension(:), intent(inout) :: wv,wmod,fmod
      real(double), dimension(:), intent(inout) :: fmodbin
   end subroutine
end interface
interface
   subroutine fitspectrum(npt,wv,flux,nmodel,wmod,fmod,fmodbin,nfit)
      use precision
      implicit none
      integer, intent(inout) :: npt,nfit,nmodel
      real(double), dimension(:), intent(inout) :: wv,flux,fmodbin,wmod,&
       fmod
   end subroutine
end interface

if(iargc().lt.2) then  !check that we have sufficient commandline arguments
   write(0,*) "Usage: specfit <spectrum.txt> <modelspec.7 [nfit]>"
   stop
endif

call getarg(1,specfile)  !get filename for photometry

nunit=10 !unit number for data spectrum
open(unit=nunit,file=specfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",specfile
   stop
endif

nmax=10000 !initital guess at the number of data points
allocate(wv(nmax),flux(nmax))
iflag=0 !flag traces data i/o
npt=0   !initalize counter for number of data points
do !we do a loop.  If there are memory errors, we can get more this way
   call readspec(nunit,nmax,npt,wv,flux,iflag) !read in spectrum
   if(iflag.eq.1) then !reallocate array space (we ran out of space)
      allocate(wv2(nmax),flux2(nmax)) !allocate temp arrays
      wv2=wv   !copy over the data we read
      flux2=flux
      deallocate(wv,flux) !deallocate data arrays
      nmax=nmax*2 !lets get more memory
      write(0,*) "warning, increasing nmax: ",nmax
      allocate(wv(nmax),flux(nmax)) !reallocate array
      do i=1,nmax/2  !copy data back into data arrays
         wv(i)=wv2(i)
         flux(i)=flux2(i)
      enddo
      deallocate(wv2,flux2) !deallocate temp arrays
      iflag=2  !set flag that we are continuing to read in data
      cycle !repeat data read loop
   endif
   exit !successively break from data read loop.
enddo
close(nunit) !close file.
!write(0,*) "npt: ",npt  !report number of data points read.

call pgopen('?') !open plotting device
call pgpage()
call pgask(.false.)
call PGPAP ( 8.0 ,1.0) !use a square 8" across
call pgsubp(1,4)
call pgpage()
call pgsch(1.5) !make the font a bit bigger
call pgslw(3)  !make the lines a bit thicker
call pgvport(0.15,0.85,0.15,0.85) !make room around the edges for labels

call plotspec(npt,wv,flux) !plot the flux.
minw=minval(wv(1:npt)) !used to make model plots on same scale
maxw=maxval(wv(1:npt))

!read in a model spectrum
!modelfile="bt-settl/lte032-5.0-0.0a+0.0.BT-Settl.spec.7"
call getarg(2,modelfile)
nunit=11 !unit number for data spectrum
open(unit=nunit,file=modelfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",modelfile
   stop
endif

nmodelmax=1500000 !initital guess at the number of data points
allocate(wmod(nmodelmax),fmod(nmodelmax))
iflag=0 !flag traces data i/o
nmodel=0   !initalize counter for number of data points
do !we do a loop.  If there are memory errors, we can get more this way
   call readmodel(nunit,nmodelmax,nmodel,wmod,fmod,iflag) !read in spectrum
   if(iflag.eq.1) then !reallocate array space (we ran out of space)
      allocate(wmod2(nmodelmax),fmod2(nmodelmax)) !allocate temp arrays
      wmod2=wmod   !copy over the data we read
      fmod2=fmod
      deallocate(wmod,fmod) !deallocate data arrays
      nmodelmax=nmodelmax*2 !lets get more memory
      write(0,*) "warning, increasing nmodelmax: ",nmodelmax
      allocate(wmod(nmodelmax),fmod(nmodelmax)) !reallocate array
      do i=1,nmodelmax/2  !copy data back into data arrays
         wmod(i)=wmod2(i)
         fmod(i)=fmod2(i)
      enddo
      deallocate(wmod2,fmod2) !deallocate temp arrays
      iflag=2  !set flag that we are continuing to read in data
      cycle !repeat data read loop
   endif
   exit !successively break from data read loop.
enddo
close(nunit) !close file.
!write(0,*) "nmodel: ",nmodel  !report number of data points read.

allocate(fmodbin(npt))
fmodbin=0.0d0
rv=0.0d0
call binmodel(npt,wv,nmodel,wmod,fmod,fmodbin,rv)
!allocate(ratio(npt))
!ratio=flux(1:npt)/fmodbin
!call pgpage()
!call plotmodel(npt,wv,ratio,minw,maxw)
!deallocate(ratio)
call pgpage()
call plotmodel(npt,wv,fmodbin,minw,maxw)

nfit=3
if(iargc().ge.3)then
   call getarg(3,cinput)
   read(cinput,*) nfit
   if(nfit.lt.1)then
      write(0,*) "nfit must be greater than 0"
      stop
   endif
endif

call fitspectrum(npt,wv,flux,nmodel,wmod,fmod,fmodbin,nfit)

call pgpage()
call plotmodel(npt,wv,fmodbin,minw,maxw)

call pgpage()
chisq=0.0d0
allocate(diff(npt))
do i=1,npt
   diff(i)=flux(i)-fmodbin(i)
   chisq=chisq+diff(i)*diff(i)
enddo
call plotmodel(npt,wv,diff,minw,maxw)
write(6,*) chisq,npt," = chi-square"

call pgclos()

end program specfit
