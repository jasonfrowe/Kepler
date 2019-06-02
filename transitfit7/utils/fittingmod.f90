module fittingmod
   use precision, only: double
   implicit none
   integer, pointer :: nbodies2f,npt2f,ntmidmax2f
   real(double), pointer :: tol2f
   real(double), dimension(:), pointer :: sol2f,time2f,flux2f,ferr2f,itime2f,y2f
   real(double), dimension(:,:), pointer :: serr2f,yserr2f
end module fittingmod
