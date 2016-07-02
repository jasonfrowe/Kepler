subroutine getfits(Refname,naxes,Ref,Rmin,Rmax,nkeys,header,bpix)
!Jason Rowe 2015 - jasonfrowe@gmail.com
use precision
implicit none
integer :: nkeys,status,unitfits,readwrite,dumi,i,nspace,nfound,npixels,&
   group,firstpix,nbuf,j,nbuffer,npt
integer, dimension(2) :: naxes
real(double) :: Rmin,Rmax,bpix,nullval
real(double), dimension(:,:) :: Ref
real(double), allocatable, dimension(:) :: buffer
character(80) :: Refname,record
character(80), dimension(:) :: header
logical :: anynull

! status will report errors.  No errors means status=0.
! initalize value of status
status=0
! gets an unused unit number to open fits file
call ftgiou(unitfits,status)
! setting to zero makes fits file readwrite
readwrite=0
! open this fits file
call ftopen(unitfits,Refname,readwrite,dumi,status)
if(status.ne.0)then
   write(0,*) "Status: ",status
   write(0,*) "Cannot open "
   write(0,'(A80)') Refname
   stop
endif

nkeys=0
! get number of headers in image
call ftghsp(unitfits,nkeys,nspace,status)
do i=1,nkeys
   call ftgrec(unitfits,i,record,status)
   header(i)=record
!   write(6,'(A80)') header(i)
enddo

call ftgknj(unitfits,'NAXIS',1,2,naxes,nfound,status)
if((naxes(1).gt.size(ref(:,1))).or.(naxes(2).gt.size(ref(1,:))))then
   write(0,*) "inadequate space for FITS."
   write(0,*) "Needed: ", naxes(1),naxes(2)
   write(0,*) "Available: ", size(ref(:,1)),size(ref(1,:))
   stop
endif
!write(0,*) naxes(1),naxes(2)

!Check that it found both NAXIS1 and NAXIS2 keywords.
if (nfound.ne.2)then
   write(6,*) 'READIMAGE failed to read the NAXISn keywords.'
   stop
endif

npixels=naxes(1)*naxes(2)
group=1
firstpix=1
nullval=bpix
nbuf=naxes(1)
j=0
allocate(buffer(nbuf))
do while (npixels.gt.0)
   nbuffer=min(nbuf,npixels)
   call ftgpvd(unitfits,group,firstpix,nbuffer,nullval,buffer,          &
    anynull,status)
   j=j+1
   do i=1,nbuffer
      Ref(i,j)=buffer(i)
   enddo
   npixels=npixels-nbuffer
   firstpix=firstpix+nbuffer
enddo

!find min-max of dataset
npt=0
do i=1,naxes(1)
   do j=1,naxes(2)
      if(Ref(i,j).lt.bpix)then !skip bad pixels
         npt=npt+1
         if(npt.eq.1)then
            Rmin=Ref(i,j)
            Rmax=Ref(i,j)
         else
            Rmin=min(Ref(i,j),Rmin)
            Rmax=max(Ref(i,j),Rmax)
         endif
      endif
   enddo
enddo


!close file
call ftclos(unitfits,status)
call ftfiou(unitfits,status)

return
end subroutine getfits
