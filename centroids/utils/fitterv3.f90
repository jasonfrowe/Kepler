subroutine fitterv3(naxes,Image,nstar,xcoo,ycoo,Ic,ngsol,gsol,sat,bpix, &
 starmap,numnei,nei,id,radnei)
use precision
use fittermodv3
implicit none
!import vars
integer :: nstar
integer, target :: ngsol
integer, dimension(2), target :: naxes
integer, dimension(:), target :: numnei
integer, dimension(:) :: id
integer, dimension(:,:), target :: starmap,nei
real(double), target :: sat
real(double) :: bpix,radnei
real(double), dimension(:) :: xcoo,ycoo,Ic,gsol
real(double), dimension(:,:), target :: Image
!local vars
integer :: i,iradneid2,iradnei,i1,i2,j1,j2,icoo,jcoo,iradneip1,    &
 npix,j,ii,nstarsub,nfitin,lwa,info
integer, target :: nfitsub,nfit
integer, dimension(2), target :: naxessub
integer, allocatable, dimension(:), target :: isolsub,numneisub,isol
integer, allocatable, dimension(:) :: iwa
integer, allocatable, dimension(:,:), target :: starmapsub,neisub
real(double) :: tol
real(double), allocatable, dimension(:), target :: solsub,sol
real(double), allocatable, dimension(:) :: solin,fvec,wa
real(double), allocatable, dimension(:,:), target :: subIm
external fcn

!precompute
iradnei=int(radnei)
iradneip1=iradnei+1
iradneid2=iradnei/2

!number of parameters to fit
nfit=ngsol+3*nstar
allocate(sol(nfit))

!initial Guesses
do i=1,ngsol
   sol(i)=gsol(i)
enddo
do i=1,nstar
   sol(ngsol+3*(i-1)+1)=Ic(i)
   sol(ngsol+3*(i-1)+2)=xcoo(i)
   sol(ngsol+3*(i-1)+3)=ycoo(i)
enddo

!allocate subarray
allocate(subIm(iradneip1,iradneip1))
!loop through all stars to improve centroids
do ii=1,nstar
   !create subraster around star of interest.
   icoo=int(xcoo(ii))
   jcoo=int(ycoo(ii))
   subIm=bpix !initialize sub-raster with bad-pixels
   i1=max(1,icoo-iradneid2)
   i2=min(naxes(1),icoo+iradneid2)
   j1=max(1,jcoo-iradneid2)
   j2=min(naxes(2),jcoo+iradneid2)
   subIm(1:i2-i1+1,1:j2-j1+1)=Image(i1:i2,j1:j2)
   naxessub(1)=i2-i1+1 !size of subraster
   naxessub(2)=j2-j1+1

   !estimate number of good pixels in subraster
   npix=0
   do i=1,iradneip1
      do j=1,iradneip1
         if(subIm(i,j).lt.sat)then
            npix=npix+1
         endif
      enddo
   enddo
!   write(0,*) "Number of pixels to fit in Image: ",npix

   !calculate number of parameters for local fit.
   nstarsub=1+numnei(ii)
   nfitsub=ngsol+3*nstarsub
   allocate(solsub(nfitsub),isolsub(nfitsub))
   do i=1,ngsol
      solsub(i)=sol(i) !shape parameters
      isolsub(i)=0     !do not fit shape parameters
   enddo
   solsub(ngsol+1)=sol(ngsol+3*(ii-1)+1) !amplitude
   isolsub(ngsol+1)=1 !fit amplitude
   solsub(ngsol+2)=sol(ngsol+3*(ii-1)+2)-dble(i1-1) !adjust x to center
   isolsub(ngsol+2)=1 !fit centroid
   solsub(ngsol+3)=sol(ngsol+3*(ii-1)+3)-dble(j1-1) !adjust y to center
   isolsub(ngsol+3)=1 !fit centroid
!   write(0,*) (solsub(i),i=1,7)

!  fill in neighbouring stars
   do i=1,numnei(ii)
      solsub(ngsol+3*(i)+1)=sol(ngsol+3*(nei(ii,i)-1)+1)
      isolsub(ngsol+3*(i)+1)=1 !fit amplitude
      solsub(ngsol+3*(i)+2)=sol(ngsol+3*(nei(ii,i)-1)+2)-dble(i1-1)
      isolsub(ngsol+3*(i)+2)=1 !fit centroid
      solsub(ngsol+3*(i)+3)=sol(ngsol+3*(nei(ii,i)-1)+3)-dble(j1-1)
      isolsub(ngsol+3*(i)+3)=1 !fit centroid
   enddo
!   write(0,*) (solsub(i),i=1,nfitsub)

   nfitin=Sum(isolsub(1:nfitsub))
   allocate(solin(nfitin))
   j=0
   do i=1,nfitsub
      if(isolsub(i).eq.1)then
         j=j+1
         solin(j)=solsub(i)
      endif
   enddo

   !allocate variables for lmdif1
   allocate(fvec(npix))
   allocate(iwa(nfitin)) !work space for fitter
   lwa=npix*nfitin+5*npix*nfitin  !work-space for fitter
   allocate(wa(lwa))

   !make a local starmap and neighbours
   allocate(starmapsub(naxessub(1),naxessub(2)),numneisub(1))
   starmapsub=1
   numneisub(1)=numnei(ii)
   allocate(neisub(1,max(1,numnei(ii))))
   do i=1,numnei(ii)
      neisub(1,i)=i+1
   enddo

   !update pointers to pass info to function for fitter
   Image2 => subIm
   sol2 => solsub
   isol2 => isolsub
   nfit2 => nfitsub
   naxes2 => naxessub
   sat2 => sat
   starmap2 => starmapsub
   ngsol2 => ngsol
   numnei2 => numneisub
   nei2 => neisub

   !parameters controling fit
   tol=1.0d-10 !fitting tolerence
   info=0     !keep track of errors from minpack

   !call fitter
!   write(0,*) "starting lmdif1"
!   write(0,*) (solin(i),i=1,3)
   call lmdif1(fcn,npix,nfitin,solin,fvec,tol,info,iwa,wa,lwa)
!   write(0,*) (solin(i),i=1,3)
!   write(0,*) "info: ",info

   j=0
   do i=1,nfitsub
      if(isolsub(i).eq.1)then
         j=j+1
         solsub(i)=solin(j)
      endif
   enddo

!   write(0,*) solsub(5),solsub(6)+dble(i1-1),solsub(7)+dble(j1-1)
!   write(0,*) Ic(ii),xcoo(ii),ycoo(ii)
   Ic(ii)=solsub(5)
   xcoo(ii)=solsub(6)+dble(i1-1)
   ycoo(ii)=solsub(7)+dble(j1-1)

!   read(5,*)

   deallocate(solsub,isolsub,fvec,iwa,wa,starmapsub,numneisub,neisub,   &
    solin)
enddo

!now we update the shape parameters.

!allocate space for isol
allocate(isol(nfit))
do i=1,ngsol
   isol(i)=1 !fit for shape parameters.
enddo

!initial Guesses
do i=1,ngsol
   sol(i)=gsol(i)
enddo
!update initial Guesses with updated centroids and peak flux.
do i=1,nstar
   sol(ngsol+3*(i-1)+1)=Ic(i) !flux
   isol(ngsol+3*(i-1)+1)=0
   sol(ngsol+3*(i-1)+2)=xcoo(i) !x-position
   isol(ngsol+3*(i-1)+2)=0
   sol(ngsol+3*(i-1)+3)=ycoo(i) !y-position
   isol(ngsol+3*(i-1)+3)=0
enddo

npix=0
do i=1,naxes(1)
   do j=1,naxes(2)
      if(Image(i,j).lt.sat)then
         npix=npix+1
      endif
   enddo
enddo
!write(0,*) "Number of pixels to fit in Image: ",npix

!number of parameters to fit
nfitin=Sum(isol(1:nfit))
!allocate space for fitter parameters
allocate(solin(nfitin))
j=0
do i=1,nfit
   if(isol(i).eq.1)then
      j=j+1
      solin(j)=sol(i)
   endif
enddo

!allocate variables for lmdif1
allocate(fvec(npix))
allocate(iwa(nfitin)) !work space for fitter
lwa=npix*nfitin+5*npix*nfitin  !work-space for fitter
allocate(wa(lwa))

!update pointers to pass info to function for fitter
Image2 => Image
sol2 => sol
isol2 => isol
nfit2 => nfit
naxes2 => naxes
sat2 => sat
starmap2 => starmap
ngsol2 => ngsol
numnei2 => numnei
nei2 => nei

!write(0,*) "Global shape fit.."
!write(0,*) solin(1:nfitin)
call lmdif1(fcn,npix,nfitin,solin,fvec,tol,info,iwa,wa,lwa)
!write(0,*) solin(1:nfitin)
!write(0,*) "info :",info

!move new fit parameters in solution.
j=0
do i=1,nfit
   if(isol(i).eq.1)then
      j=j+1
      sol(i)=solin(j)
   endif
enddo

!initial Guesses
do i=1,ngsol
   gsol(i)=sol(i)
enddo

return
end subroutine fitterv3

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn(m,n,x,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittermodv3
implicit none
integer :: m,n,iflag
real(double), dimension(n) :: x
real(double), dimension(m) :: fvec
!local vars
integer :: i,j,k,ii,npix,kk
real(double) :: di,dj,gpixmod,chisq
real(double), dimension(7) :: lsol
real(double), allocatable, dimension(:) :: sol

allocate(sol(nfit2))
sol=sol2
j=0
do i=1,nfit2
   if(isol2(i).eq.1)then
      j=j+1
      sol(i)=x(j)
!      write(0,*) i,j,sol2(i),sol(i)
   endif
enddo
!write(0,*) "sol: ",(sol(i),i=1,nfit2)
!read(5,*)

!n is the number of fitted variables
lsol(1)=sol(1)
lsol(2)=sol(2)
lsol(3)=sol(3)
lsol(4)=sol(4)

!write(0,*) "lsol: ",(lsol(i),i=1,4)

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
         lsol(5)=abs(sol(ngsol2+3*(k-1)+1)) !amplitude
         lsol(6)=sol(ngsol2+3*(k-1)+2)-di !dx
         lsol(7)=sol(ngsol2+3*(k-1)+3)-dj !dy
!         write(0,*) "lsol: ",(lsol(ii),ii=1,7)
         fvec(npix)=fvec(npix)+gpixmod(lsol)
!         write(0,*) i,j,fvec(npix),Image2(i,j)
!         read(5,*)
         do ii=1,numnei2(k) !deal with neighbours
            kk=nei2(k,ii)
            lsol(5)=abs(sol(ngsol2+3*(kk-1)+1)) !amplitude
            lsol(6)=sol(ngsol2+3*(kk-1)+2)-di !dx
            lsol(7)=sol(ngsol2+3*(kk-1)+3)-dj !dy
            fvec(npix)=fvec(npix)+gpixmod(lsol)
!            write(0,*) "gsol: ",gpixmod(lsol)
         enddo
!         write(0,*) "lsol: ",(lsol(ii),ii=1,7)
!         write(0,*) i,j,fvec(npix),Image2(i,j)
!         read(5,*)
         fvec(npix)=(fvec(npix)-Image2(i,j))/sqrt(abs(Image2(i,j)))
      endif
   enddo
enddo

!chisq=0.0d0
!do i=1,npix
!   chisq=chisq+fvec(i)*fvec(i)
!enddo
!write(0,*) "chisq: ",chisq
!read(5,*)

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
