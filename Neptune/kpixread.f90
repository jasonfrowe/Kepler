program kpixread
use precision
implicit none
integer :: iargc,nunit,filestatus,nfile,nsteps,nstep,nfilemax,i,npixmax, &
   npix,j,nxmax,nymax,k,nkeys,nkeysmax,nunitc,flag,nskip,iter,ix,iy,ii
integer, dimension(2) :: crval
real(double) :: tpix,bpix,tavg,tavg1,x,y
real(double), allocatable, dimension(:) :: pix,xcoo,ycoo
real(double), allocatable, dimension(:,:) :: parray,tarray,refarray
character(80) :: listfile,fitsfile,gfile
character(80), allocatable, dimension(:) :: filenames, header

!below are all the interfaces to allow dynamic arrays.
interface
   subroutine fitsread(fitsfile,nstep,npix,tpix,pix,crval,nkeysmax,     &
    nkeys,header,ix,iy)
      use precision
      implicit none
      integer, intent(inout) :: nstep,npix,nkeysmax,nkeys,ix,iy
      integer, dimension(2) :: crval
      real(double) :: tpix
      real(double), dimension(:), intent(inout) :: pix
      character(80), intent(inout) :: fitsfile
      character(80), dimension(:), intent(inout) :: header
   end subroutine fitsread
end interface
interface
   subroutine getnsteps(fitsfile,nsteps,npixmax)
      use precision
      implicit none
      integer, intent(inout) :: nsteps,npixmax
      character(80), intent(inout) :: fitsfile
   end subroutine getnsteps
end interface
interface
   subroutine displayfits(nxmax,nymax,parray,bpix,tavg,refarray)
      use precision
      implicit none
      integer, intent(inout) :: nxmax,nymax
      real(double), dimension(:,:), intent(inout) :: parray,refarray
      real(double), intent(inout) :: bpix,tavg
   end subroutine displayfits
end interface
interface
   subroutine writefits(nxmax,nymax,parray,bpix,tavg,nkeys,header,nstep)
      use precision
      implicit none
      integer :: nxmax,nymax,nkeys,nstep
      real(double) :: bpix,tavg
      real(double), dimension(:,:) :: parray
      character(80), dimension(:) :: header
   end subroutine writefits
end interface

open(unit=11,file='minmax.dat')

if(iargc().lt.1)then
   write(0,*) "Usage: kpixread input.list"
   stop
endif

call getarg(1,listfile)  !get filename for photometry

nunit=10 !unit number for file list
open(unit=nunit,file=listfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",listfile
   stop
endif

!read in filenames
nfilemax=1000
allocate(filenames(nfilemax))
nfile=0
do
   read(nunit,*,iostat=filestatus) fitsfile
   if(filestatus == 0) then
      nfile=nfile+1 !count number of files
      if(nfile.eq.1)then
         call getnsteps(fitsfile,nsteps,npixmax) !predetermine number of time-steps
         write(0,*) "NSteps: ",nsteps,npixmax
!         read(5,*)
      endif
      filenames(nfile)=fitsfile
      cycle
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!!",listfile
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit)
write(0,*) "nfiles: ",nfile

allocate(xcoo(nsteps),ycoo(nsteps))
xcoo=10.0 !make initial centroids 'bad'
ycoo=10.0
!nunitc=12 !unit number for file list
!open(unit=nunitc,file="centroids.dat",iostat=filestatus,status='old')
!if(filestatus>0)then !trap missing file errors
!   write(0,*) "Cannot open centroids.dat"
!   stop
!endif
!do
!   read(nunitc,*,iostat=filestatus) nstep,x,y
!   if(filestatus == 0) then
!      xcoo(nstep)=x
!      ycoo(nstep)=y
!      cycle
!   elseif(filestatus == -1) then
!      exit  !successively break from data read loop.
!   else
!      write(0,*) "File Error!! centroids.dat"
!      write(0,'(A8,I3)') "iostat: ",filestatus
!      stop
!   endif
!enddo
!close(nunitc)

!call pgopen('/xwindow')
!cp call pgopen('movie.ps/vcps')
!call PGPAP ( 8.0 ,9.0/16.0) !use a square 8" across
!call pgpage()

nxmax=1024*4
nymax=1024*4
nkeysmax=700
allocate(tarray(nxmax,nymax),parray(nxmax,nymax),header(nkeysmax),      &
   refarray(nxmax,nymax))
bpix=1000000.00 !bad pixel value

!wasn't npixmax pre-determined?
npixmax=10000 ! buffer for reading in pixels
allocate(pix(npixmax))
!loop over all files in listfile
nstep=1
!nstep=30000
iter=1
do while(nstep.lt.nsteps)
!do while(nstep.lt.82520)
!   write(0,*) "xcoo,ycoo: ",xcoo(nstep),ycoo(nstep)
   flag=0
   nskip=0
!   do while(flag.eq.0)
!      if((xcoo(nstep).lt.0.3).and.(xcoo(nstep).gt.-0.2).and.            &
!       (ycoo(nstep).lt.0.3).and.(ycoo(nstep).gt.-0.2))then
!         flag=1
!      else
!         nstep=nstep+1
!         nskip=nskip+1
!      endif
!   enddo
   if(nstep.lt.10)then
      write(gfile,501) 'p00000',nstep,'.png/png'
      501 format(A6,I1,A8)
   elseif(nstep.lt.100)then
      write(gfile,502) 'p0000',nstep,'.png/png'
      502 format(A5,I2,A8)
   elseif(nstep.lt.1000)then
      write(gfile,503) 'p000',nstep,'.png/png'
      503 format(A4,I3,A8)
   elseif(nstep.lt.10000)then
      write(gfile,504) 'p00',nstep,'.png/png'
      504 format(A3,I4,A8)
   elseif(nstep.lt.100000)then
      write(gfile,505) 'p0',nstep,'.png/png'
      505 format(A2,I5,A8)
   else
      write(gfile,506) 'p',nstep,'.png/png'
      506 format(A1,I6,A8)
   endif
!   write(0,*) gfile
   call pgopen(gfile)
   call PGPAP ( 22.5761 ,9.0/16.0) !use a square 8" across
   call pgpage()

   parray=bpix+1.0
   tarray=0.0d0
   tavg=0.0d0
   do i=1,nfile
      fitsfile=filenames(i)
!      write(0,500) fitsfile
      500 format(A80)
      call fitsread(fitsfile,nstep,npix,tpix,pix,crval,nkeysmax,nkeys,  &
         header,ix,iy)
!      write(0,*) crval(1),crval(2),ix,iy
!      read(5,*)
!      write(0,*) tpix,crval(1),crval(2)
!      write(0,*) (pix(j),j=1,npix)

!      do j=1,npix
!         tarray(crval(1),crval(2)+j-1)=tpix
!         parray(crval(1),crval(2)+j-1)=pix(j)
!      enddo
      ii=0
      do k=1,iy
         do j=1,ix
            ii=ii+1
            parray(crval(1)+j-1,crval(2)+k-1)=pix(ii)
         enddo
      enddo
!      write(0,*) i,tpix
      tavg=tavg+tpix
!      read(5,*)
   enddo
   tavg=tavg/dble(nfile)
   if(iter.eq.1) then
      tavg1=tavg
      write(0,*) "Tavg1: ",tavg1
   endif
   tavg=tavg-tavg1

!   do i=1,nkeys
!      write(0,'(A80)') header(i)
!   enddo

!   write(0,*) gfile

!  transpose Image for movie
!   do i=1,nxmax
!      do j=1,nymax
!         tarray(j,i)=parray(i,j)
!      enddo
!   enddo
!   parray=tarray

   if(iter.eq.1)then
      write(0,*) "copying reference"
      refarray=parray
   endif
   call displayfits(nxmax,nymax,parray,bpix,tavg,refarray)

   call writefits(nxmax,nymax,parray,bpix,tavg,nkeys,header,nstep)
!   read(5,*)

!   if(nstep.eq.1)then
!      nstep=24081
!   else
      nstep=nstep+1-nskip
!   endif
   iter=iter+1
   call pgclos()

enddo

close(11)

!call pgclos()

end program kpixread
