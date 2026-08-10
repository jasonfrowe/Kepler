subroutine cutoutliers(npt,x,y,yerr,itime)
use precision
implicit none
!import variable
integer :: npt
real(double), dimension(:) :: x,y,yerr,itime
!local variables
integer :: i,j,i1,i2,nsamp,nsampmax
integer, allocatable, dimension(:) :: icut
real(double) :: threshold,vp,vm,sigma,std,meddiff
real(double), allocatable, dimension(:) :: samps

!old thresholding variable, this value is overwriten if 'sigma > 0'
threshold=0.0005 !Neptune

!new thresholding vars.  Only works if sigma > 0
nsampmax=3 !number of +/- nearby samples to use for stats
sigma=3.0

allocate(icut(npt),samps(nsampmax*2+1))
icut=0 !intialize to keep all data by default

!$OMP PARALLEL DO PRIVATE(vp,vm,i1,i2,nsamp,samps,std,threshold)
do i=2,npt-1

   !estimate threshold by taking a local sample, then estimating the standard deviation  
   i1=max(1,i-nsampmax)
   i2=min(npt,i+nsampmax)
   nsamp=i2-i1+1
   samps(1:nsamp)=y(i1:i2)
   std=meddiff(nsamp,samps) !calculate the median point-to-point change.  
   if(sigma.gt.0)threshold=std*sigma

   vp=y(i)-y(i+1)
   vm=y(i)-y(i-1)
   if((abs(vp).gt.threshold).and.(abs(vm).gt.threshold).and.(vp/vm.gt.0))then
      icut(i)=1 !cut
   endif

   !if(x(i)+2.5653031294d1.gt.28.1)then
   !   write(0,'(4(F12.6,1X),I2)') x(i)+2.5653031294d1,std,vp,vm,icut(i)
   !   read(5,*)
   !endif

enddo
!$OMP END PARALLEL DO

j=0
do i=1,npt
   if(icut(i).eq.0)then
      j=j+1
      x(j)=x(i)
      y(j)=y(i)
      yerr(j)=yerr(i)
      itime(j)=itime(i)
   endif
enddo
npt=j

return
end
