subroutine fitter(naxes,Image,nstar,Ic,xcoo,ycoo,ngsol,gsol,sat,bpix,   &
 starmap,numnei,nei)
use precision
use fittermod
implicit none
!import vars
integer :: nstar
integer, target :: ngsol
integer, dimension(2), target :: naxes
integer, dimension(:), target :: numnei
integer, dimension(:,:), target :: starmap,nei
real(double) :: bpix
real(double), target :: sat
real(double), dimension(:) :: Ic,xcoo,ycoo,gsol
real(double), dimension(:,:), target :: Image
!local vars
integer i,j,npix,lwa,nfit,info
integer, allocatable, dimension(:) :: iwa
real(double) :: tol
real(double), allocatable, dimension(:) :: sol,fvec,wa
external fcn

npix=0
do i=1,naxes(1)
   do j=1,naxes(2)
      if(Image(i,j).lt.sat)then
         npix=npix+1
      endif
   enddo
enddo
write(0,*) "Number of pixels to fit in Image: ",npix

!number of parameters to fit
nfit=ngsol+3*nstar
allocate(sol(nfit))

!initial Guesses
sol(1)=0.0
sol(2)=1.3
sol(3)=1.3
sol(4)=0.0
do i=1,nstar
   sol(ngsol+3*(i-1)+1)=Image(int(xcoo(i)),int(ycoo(i)))*0.6
   sol(ngsol+3*(i-1)+2)=xcoo(i)
   sol(ngsol+3*(i-1)+3)=ycoo(i)
enddo

!allocate variables for lmdif1
allocate(fvec(npix))
allocate(iwa(nfit)) !work space for fitter
lwa=npix*nfit+5*npix*nfit  !work-space for fitter
allocate(wa(lwa))

!update pointer to pass info to function for fitter
Image2 => Image
naxes2 => naxes
sat2 => sat
starmap2 => starmap
nei2 => nei
numnei2 => numnei
ngsol2 => ngsol

!parameters controling fit
tol=1.0d-10 !fitting tolerence
info=0     !keep track of errors from minpack

!call fitter
write(0,*) "starting lmdif1"
call lmdif1(fcn,npix,nfit,sol,fvec,tol,info,iwa,wa,lwa)
write(0,*) "info: ",info
write(0,*) (sol(i),i=1,10)

return
end subroutine fitter


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn(m,n,x,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittermod
implicit none
integer :: m,n,iflag
real(double), dimension(n) :: x
real(double), dimension(m) :: fvec
!local vars
integer :: i,j,k,ii,npix,kk
real(double) :: di,dj,gpixmod,chisq
real(double), dimension(7) :: lsol

!n is the number of fitted variables
lsol(1)=x(1)
lsol(2)=x(2)
lsol(3)=x(3)
lsol(4)=x(4)

!write(0,*) "x: ",(x(i),i=1,7)

fvec=0.0d0 !Initialize to zero.

npix=0
do i=1,naxes2(1)
   di=dble(i)
   do j=1,naxes2(2)
      dj=dble(j)
      if(Image2(i,j).lt.sat2)then
         npix=npix+1
         if(npix.gt.m)then
            write(0,*) "Npix greater than m.. "
            write(0,*) "number of pixels not counted correctly"
            stop
         endif
         k=starmap2(i,j)
         lsol(5)=x(ngsol2+3*(k-1)+1) !amplitude
         lsol(6)=x(ngsol2+3*(k-1)+2)-di !dx
         lsol(7)=x(ngsol2+3*(k-1)+3)-dj !dy
         fvec(npix)=fvec(npix)+gpixmod(lsol)
!         write(0,*) (lsol(ii),ii=1,7)
!         write(0,*) fvec(npix),Image2(i,j)
!         read(5,*)
         do ii=1,numnei2(k) !deal with neighbours
            kk=nei2(k,ii)
            lsol(5)=x(ngsol2+3*(kk-1)+1) !amplitude
            lsol(6)=x(ngsol2+3*(kk-1)+2)-di !dx
            lsol(7)=x(ngsol2+3*(kk-1)+3)-dj !dy
            fvec(npix)=fvec(npix)+gpixmod(lsol)
         enddo
         fvec(npix)=(fvec(npix)-Image2(i,j))/sqrt(abs(Image2(i,j)))
      endif
   enddo
enddo

chisq=0.0d0
do i=1,npix
   chisq=chisq+fvec(i)*fvec(i)
enddo
write(0,*) "chisq: ",chisq

return
end subroutine fcn

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function gpixmod(sol)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
real(double) :: gpixmod
real(double), dimension(3) :: temp
real(double), dimension(7) :: sol

temp(1)=-sol(6)*sol(6)/(sol(2)*sol(2))
temp(2)=-sol(7)*sol(7)/(sol(3)*sol(3))
temp(3)=2.0d0*sol(4)*sol(6)*sol(7)/(sol(2)*sol(3))

gpixmod=sol(1)+sol(5)*exp(temp(1)+temp(2)+temp(3))

return
end function gpixmod

