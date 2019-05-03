module fittingmod
   use precision, only: double
   implicit none
   integer isol,im,iy
   integer, pointer :: nbodies2,npt2
   real(double), pointer :: tol2
   real(double), dimension(:), pointer :: sol2,y2,time2,flux2,ferr2,itime2
   real(double), dimension(:,:), pointer :: serr2,yserr2
end module fittingmod
