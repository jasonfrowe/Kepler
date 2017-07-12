subroutine binmodel(npt,wv,nmodel,wmod,fmod,fmodbin,rv)
use precision
implicit none
integer :: npt,nmodel,i,nbin
integer, allocatable, dimension(:) :: q,nfmodbin
real(double), parameter :: c=2.99792458d8
real(double) :: wvdiffmed,zpt,rv,dl
real(double), dimension(:) :: wv,wmod,fmod,fmodbin
real(double), allocatable, dimension(:) :: wvdiff

!calculate spacing between observed spectra.. assumes spectra is sorted short to long
allocate(q(npt-1),wvdiff(npt-1)) !allocate work arrays
do i=1,npt-1
   wvdiff(i)=wv(i+1)-wv(i) !calculate difference
enddo
call rqsort(npt-1,wvdiff,q) !sort differences
wvdiffmed=wvdiff(q((npt-1)/2)) !get median value
!write(0,*) "wvdiffmed: ",wvdiffmed
deallocate(q,wvdiff) !do not need work arrays
!write(0,*) "deallocated q,wvdiff"
zpt=wv(1)-wvdiffmed/2.0d0
!write(0,*) "zpt: ",zpt

!write(0,*) "allocating"
allocate(nfmodbin(npt))
!write(0,*) "done allocating"
nfmodbin=0
 fmodbin=0.0d0
do i=1,nmodel
!   write(0,*) i
   dl=wmod(i)*rv/c
!   write(0,*) dl
   nbin=int((wmod(i)-zpt+dl)/wvdiffmed)+1
!   write(0,*) dl,nbin
   if((nbin.gt.0).and.(nbin.le.npt))then
      nfmodbin(nbin)=nfmodbin(nbin)+1
       fmodbin(nbin)= fmodbin(nbin)+fmod(i)
   endif
enddo
do i=1,npt
   fmodbin(i)=fmodbin(i)/dble(nfmodbin(i))
enddo
!deallocate(nfmodbin)

return
end subroutine binmodel
