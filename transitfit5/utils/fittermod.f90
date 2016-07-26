module fittermod
   use precision, only: double
   implicit none
   integer, pointer:: nfit2,nplanet2,nfrho2
   integer, dimension(:), pointer :: dtype2,ntt2
   real(double), pointer :: rhoi2
   real(double), dimension(:), pointer :: sol2,aT2,aM2,aE2,aIT2,rhoierr2
   real(double), dimension(:,:), pointer :: serr2,tobs2,omc2
end module fittermod
