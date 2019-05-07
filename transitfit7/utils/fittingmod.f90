module fittingmod
   use precision, only: double
   implicit none
   integer isol,im,iy
   integer, pointer :: nbodies2,npt2,ntmidmax2
   integer, dimension(:), pointer :: ntmid2
   real(double), pointer :: tol2
   real(double), dimension(:), pointer :: sol2,time2,flux2,ferr2,itime2,y2
   real(double), dimension(:,:), pointer :: serr2,yserr2,tmid2
end module fittingmod
