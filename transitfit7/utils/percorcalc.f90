subroutine percorcalc(nbodies,sol,ntmidmax,ntmid,tmid,percor)
use precision
implicit none
!import vars
integer :: nbodies,ntmidmax
integer, dimension(:) :: ntmid
real(double), dimension(:) :: sol,percor
real(double), dimension(:,:) :: tmid
!local vars
integer :: i,j,k,np
real(double) :: per,ptest
!real(double) :: median
real(double), allocatable, dimension(:) :: p
 
allocate(p(ntmidmax)) !stored period from n-body positions

percor=0.0d0 !default is to not have a period correction 

!loop for all nbodies, except central source
do i=2,nbodies
  np=7+7*(i-1)
  Per=sol(np+2) !model period 
  !scan through all the measured timings
  k=0 !counts number of complete nbody-periods found
  if (ntmid(i).gt.2) then  !need at least three measurements to get a mean period
  	do j=2,ntmid(i)
  		ptest = tmid(i,j) - tmid(i,j-1) !the nbody period
  		!the nbody period should be close to the wanted period...
  		if (abs((ptest-Per)/Per).lt.0.50) then  !if greater than 1, then a transit was missed. 
  			k=k+1
  			p(k)=ptest
       !write(0,502) k,tmid(i,j),ptest,Per
        502 format(I3,1X,28(1PE15.8,1X))
  		endif 
  	enddo
    if (k.gt.1) then !make sure we have periods..
      percor(i)=Sum(p(1:k))/dble(k)-Per  !mean
     !percor(i)=median(k,p)-Per
    endif
  endif
  !write(0,501) "percor",i,ntmid(i),percor(i),Per,percor(i)+Per
  !read(5,*)
  501 format(A6,1X,I3,1X,I3,1X,28(1PE15.8,1X))
enddo



return
end
