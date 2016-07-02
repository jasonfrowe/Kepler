program psffit
use precision
implicit none
integer :: iargc,xmax,ymax,nkeys,nkeysmax,nstarmax,nstar,ngsol,neimax,i,&
 niter,j
integer, dimension(2) :: naxes
integer, allocatable, dimension(:) :: id,numnei
integer, allocatable, dimension(:,:) :: starmap,nei
real(double) :: Rmin,Rmax,bpix,sat,radnei
real(double), allocatable, dimension(:) :: xcoo,ycoo,Ic,gsol
real(double), allocatable, dimension(:,:) :: Image,tImage
character(80) :: fitsfile,coofile
character(80), allocatable, dimension(:) :: header

!Interfaces to subroutines
interface
   subroutine getfits(Refname,naxes,Ref,Rmin,Rmax,nkeys,header,bpix)
      use precision
      implicit none
      integer :: nkeys
      integer, dimension(2), intent(inout) :: naxes
      real(double), intent(inout) :: Rmin,Rmax,bpix
      real(double), dimension(:,:), intent(inout) :: Ref
      character(80), intent(inout) :: Refname
      character(80), dimension(:), intent(inout) :: header
   end subroutine getfits
end interface
interface
   subroutine fitterv3(naxes,Image,nstar,xcoo,ycoo,Ic,ngsol,gsol,sat,   &
    bpix,starmap,numnei,nei,id,radnei)
      use precision
      implicit none
      integer :: nstar
      integer, target :: ngsol
      integer, dimension(2), target :: naxes
      integer, dimension(:), target :: numnei
      integer, dimension(:) :: id
      integer, dimension(:,:), target :: starmap,nei
      real(double) :: bpix,radnei
      real(double), target :: sat
      real(double), dimension(:) :: xcoo,ycoo,Ic,gsol
      real(double), dimension(:,:), target :: Image
   end subroutine fitterv3
end interface
interface
   subroutine cstarmap(nstar,xcoo,ycoo,naxes,Image,starmap,bpix)
      use precision
      implicit none
      integer :: nstar
      integer, dimension(2) :: naxes
      integer, dimension(:,:) :: starmap
      real(double) :: bpix
      real(double), dimension(:) :: xcoo, ycoo
      real(double), dimension(:,:) :: Image
   end subroutine cstarmap
end interface
interface
   subroutine findnei(nstar,xcoo,ycoo,neimax,numnei,nei,radnei)
      use precision
      implicit none
      integer :: nstar,neimax
      integer, dimension(:) :: numnei
      integer, dimension(:,:) :: nei
      real(double) :: radnei
      real(double), dimension(:) :: xcoo,ycoo
   end subroutine findnei
end interface
interface
   subroutine fluxcoos(nstar,xcoo,ycoo,naxes,Image,bpix)
      use precision
      implicit none
      integer :: nstar
      integer, dimension(2) :: naxes
      real(double) :: bpix
      real(double), dimension(:) :: xcoo,ycoo
      real(double), dimension(:,:) :: Image
   end subroutine fluxcoos
end interface

!options
bpix=1000000.0 !badpixel value
sat=100000.0 !satuation threshold
xmax=2048 !maximum dimension for image (x)
ymax=2048 !maximum dimension for image (y)
nkeysmax=700 !number of lines in header we can use
nstarmax=10000 !maximum stars that can be handled
neimax=20 !maximum number of stellar neighbours
radnei=30.0 !radius to find stellar neighbours
niter=10 !number of interactions to improve PSF fits
ngsol=4 !number of global parameters for PSF shape
allocate(gsol(ngsol))
gsol(1) = 0.0 !Ib - background
gsol(2) = 1.3 !Sx - x-width
gsol(3) = 1.3 !Sy - y-width
gsol(4) = 0.0 !Sxy - xy-width

!check commandline arguments
if(iargc().lt.2)then
   write(0,*) "Usage: psffit <FITS> <COOFILE>"
   write(0,*) " <FITS> : FITS image"
   write(0,*) " <COOFILE> : file containing ID,X,Y positions"
   stop
endif

!readin commandline arguments
call getarg(1,fitsfile)
call getarg(2,coofile)

!allocate space for Image
allocate(Image(xmax,ymax))
!allocate space for header
allocate(header(nkeysmax))

!read in FITS file
call getfits(fitsfile,naxes,Image,Rmin,Rmax,nkeys,header,bpix)
write(0,*) "Rmin,Rmax: ",Rmin,Rmax
write(0,*) "naxes ",naxes(1),naxes(2)
!compact array
allocate(tImage(xmax,ymax))
tImage=Image
deallocate(Image)
xmax=naxes(1)
ymax=naxes(2)
allocate(Image(xmax,ymax))
Image(1:xmax,1:ymax)=tImage(1:xmax,1:ymax)
deallocate(tImage)

!allocate space to read in star info
allocate(id(nstarmax),xcoo(nstarmax),ycoo(nstarmax),Ic(nstarmax))

!read in co-ordinates of stars
call readstars(coofile,nstarmax,nstar,id,xcoo,ycoo)
write(0,*) "Stars read: ",nstar

!update positions with flux weighted centroids
call fluxcoos(nstar,xcoo,ycoo,naxes,Image,bpix)
!get initial estimate of Ic
do i=1,nstar
   Ic(i)=Image(int(xcoo(i)),int(ycoo(i)))*0.6
enddo


!make map of closest star for each pixel
allocate(starmap(xmax,ymax))
call cstarmap(nstar,xcoo,ycoo,naxes,Image,starmap,bpix)

!find stellar neighbours
allocate(nei(nstar,neimax),numnei(nstar))
call findnei(nstar,xcoo,ycoo,neimax,numnei,nei,radnei)

write(0,*) "Starting fitter.."
do j=1,niter
   !call datafitter
   call fitterv3(naxes,Image,nstar,xcoo,ycoo,Ic,ngsol,gsol,sat,bpix,    &
    starmap,numnei,nei,id,radnei)
!   write(0,500) (gsol(i),i=1,ngsol)
   write(0,500) Ic(1),xcoo(1),ycoo(1)
   500 format(4(1PE17.10,1X))
enddo
write(0,*) "Done fitting.."

call exportsol(nstar,Id,Ic,xcoo,ycoo,ngsol,gsol)

end program psffit
