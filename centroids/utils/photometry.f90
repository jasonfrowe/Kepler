subroutine photometry(filelist,nstar,xcoo,ycoo,naxes,napmap,apmap,&
 phot,nfilemax,bpix,nkeysmax,time,np)
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
!local vars
integer :: nunit,filestatus,i,j,k,nkeys,ii
integer, dimension(2) :: naxesin
real(double) :: Rmin,Rmax,t1
real(double), allocatable, dimension(:,:) :: Image
character(80) :: fitsfile
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

nunit=10
open(unit=nunit,file=filelist,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filelist
   stop
endif

!allocate space for image
allocate(Image(naxes(1),naxes(2)))
!allocate space for header
allocate(header(nkeysmax))

!initialize photometry array to zero
phot=0.0d0

i=0
do
   if(i.gt.nfilemax)then
      write(0,*) "Increase nfilemax to match files to read in"
      write(0,*) "nfilemax: ",nfilemax
      stop
   endif
   read(nunit,*,iostat=filestatus) fitsfile
   if(filestatus == 0) then
      i=i+1
!      write(0,*) "filename: ",fitsfile
      !read in the FITS data
      call getfits(fitsfile,naxesin,Image,Rmin,Rmax,nkeys,header,bpix)
      call gettime(t1,nkeys,header)
      time(i)=t1 !same time-stamp
      do j=1,naxes(1)
         do k=1,naxes(2)
            if((napmap(j,k).gt.0).and.(Image(j,k).lt.bpix))then
               do ii=1,napmap(j,k) !photometry!
                  phot(apmap(j,k,ii),i)=phot(apmap(j,k,ii),i)+Image(j,k)
               enddo
            endif
         enddo
      enddo
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
np=i !number of photometric measurements per star

return
end subroutine photometry

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine gettime(t1,nkeys,header)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nkeys
real(double) :: t1
character(80), dimension(nkeys) :: header
!local vars
integer :: i

do i=1,nkeys
   if(header(i)(1:8).eq."TIMEOBS")then
      read(header(i)(10:26),*) t1
!      write(0,*) "t1: ",t1
   endif
enddo

return
end subroutine gettime
