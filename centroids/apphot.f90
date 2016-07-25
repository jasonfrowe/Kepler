program makeapmask
use precision
implicit none
integer :: iargc,xmax,ymax,nkeysmax,nstarmax,nkeys,nstar,neimax,        &
 nfilemax,np
integer, dimension(2) :: naxes
integer, allocatable, dimension(:) :: id
integer, allocatable, dimension(:,:) :: napmap
integer, allocatable, dimension(:,:,:) :: apmap
real(double) :: sat,bpix,Rmin,Rmax
real(double), allocatable, dimension(:) :: xcoo,ycoo,time
real(double), allocatable, dimension(:,:) :: Image,tImage,phot
character(80) :: fitsfile,coofile,filelist
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
interface
   subroutine makeapmap(nstar,xcoo,ycoo,naxes,Image,bpix,neimax,napmap, &
    apmap)
      use precision
      implicit none
      integer :: nstar,neimax
      integer, dimension(2) :: naxes
      integer, dimension(:,:) :: napmap
      integer, dimension(:,:,:) :: apmap
      real(double) :: bpix
      real(double), dimension(:) :: xcoo,ycoo
      real(double), dimension(:,:) :: Image
   end subroutine makeapmap
end interface
interface
   subroutine photometry(filelist,nstar,xcoo,ycoo,naxes,napmap,   &
    apmap,phot,nfilemax,bpix,nkeysmax,time,np)
      use precision
      implicit none
      integer :: nstar,nfilemax,nkeysmax,np
      integer, dimension(2) :: naxes
      integer, allocatable, dimension(:,:) :: napmap
      integer, allocatable, dimension(:,:,:) :: apmap
      real(double) :: bpix
      real(double), dimension(:) :: xcoo,ycoo,time
      real(double), dimension(:,:) :: phot
      character(80) :: filelist
   end subroutine photometry
end interface
interface
   subroutine exportphot(np,nstar,phot,id,time)
      use precision
      implicit none
      integer :: np,nstar
      integer, dimension(:) :: id
      real(double), dimension(:) :: time
      real(double), dimension(:,:) :: phot
   end subroutine exportphot
end interface


!options
bpix=1000000.0 !badpixel value
sat=100000.0 !satuation threshold
xmax=2048 !maximum dimension for image (x)
ymax=2048 !maximum dimension for image (y)
nkeysmax=700 !number of lines in header we can use
nstarmax=10000 !maximum stars that can be handled
neimax=20 !maximum number of stellar neighbours
nfilemax=2000 !maximum number of files

!check commandline arguments
if(iargc().lt.3)then
   write(0,*) "Usage: psffit <FITS> <COOFILE> <LIST>"
   write(0,*) " <FITS> : reference FITS image to build ap mask"
   write(0,*) " <COOFILE> : file containing ID,X,Y positions"
   write(0,*) " <LIST>: file containing names of files for photometry"
   stop
endif

!readin commandline arguments
call getarg(1,fitsfile)
call getarg(2,coofile)
call getarg(3,filelist)

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
allocate(id(nstarmax),xcoo(nstarmax),ycoo(nstarmax))

!read in co-ordinates of stars
call readstars(coofile,nstarmax,nstar,id,xcoo,ycoo)
write(0,*) "Stars read: ",nstar

!update positions with flux weighted centroids
call fluxcoos(nstar,xcoo,ycoo,naxes,Image,bpix)

!allocate space for ap map
allocate(apmap(xmax,ymax,neimax),napmap(xmax,ymax))
!create ap masks for each star
call makeapmap(nstar,xcoo,ycoo,naxes,Image,bpix,neimax,napmap,apmap)

!allocate space
allocate(phot(nstar,nfilemax),time(nfilemax))

call photometry(filelist,nstar,xcoo,ycoo,naxes,napmap,apmap,phot,       &
 nfilemax,bpix,nkeysmax,time,np)

call exportphot(np,nstar,phot,id,time)

end program makeapmask

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine exportphot(np,nstar,phot,id,time)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: np,nstar
integer, dimension(:) :: id
real(double), dimension(:) :: time
real(double), dimension(:,:) :: phot
!local vars
integer :: i,j,nunit
character(80) :: filename

nunit=10

do i=1,nstar
!   if(id(i).eq.10002261)then
   if(id(i).gt.99999999)then
      write(filename,'(A3,I9,A4)') "klc",id(i),".dat"
   elseif(id(i).gt.9999999)then
      write(filename,'(A4,I8,A4)') "klc0",id(i),".dat"
   elseif(id(i).gt.999999)then
      write(filename,'(A5,I7,A4)') "klc00",id(i),".dat"
   elseif(id(i).gt.99999)then
      write(filename,'(A6,I6,A4)') "klc000",id(i),".dat"
   elseif(id(i).gt.9999)then
      write(filename,'(A7,I5,A4)') "klc0000",id(i),".dat"
   elseif(id(i).gt.999)then
      write(filename,'(A8,I4,A4)') "klc00000",id(i),".dat"
   elseif(id(i).gt.99)then
      write(filename,'(A9,I3,A4)') "klc000000",id(i),".dat"
   elseif(id(i).gt.9)then
      write(filename,'(A10,I2,A4)') "klc0000000",id(i),".dat"
   else
      write(filename,'(A11,I1,A4)') "klc00000000",id(i),".dat"
   endif
   open(unit=nunit,file=filename)
   do j=1,np
      if(phot(i,j).gt.0.0) write(nunit,500) time(j),phot(i,j),sqrt(phot(i,j)),j
   enddo
   close(nunit)
enddo

500 format(3(1PE17.10,1X),I6)

return
end subroutine exportphot
