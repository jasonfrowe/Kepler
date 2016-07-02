module fittermodv3
   use precision, only: double
   implicit none
   integer, pointer :: nfit2,ngsol2
   integer, dimension(:), pointer :: isol2,naxes2,numnei2
   integer, dimension(:,:), pointer :: starmap2,nei2
   real(double), pointer :: sat2
   real(double), dimension(:), pointer :: sol2
   real(double), dimension(:,:), pointer :: Image2
end module fittermodv3
