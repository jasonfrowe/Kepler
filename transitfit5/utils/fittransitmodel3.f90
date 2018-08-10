subroutine fittransitmodel3(nfit,sol,serr,nplanet,npt,aT,aM,aE,aIT,     &
 dtype,ntt,tobs,omc,nfrho,rhoi,rhoierr)
use precision
use fittermod
implicit none
integer, target :: nfit,npt,nplanet,nfrho
integer, dimension(:), target :: dtype,ntt
real(double), target :: rhoi
real(double), dimension(:), target :: sol,aT,aM,aE,aIT,rhoierr
real(double), dimension(:,:), target :: serr,tobs,omc
!local vars
integer :: i
real(double), allocatable, dimension(:) :: solin
!lmdif vars
integer :: lwa,m,n,info
real(double) :: tol
integer, allocatable, dimension(:) :: iwa
real(double), allocatable, dimension(:) :: wa,fvec
external fcn

allocate(solin(nfit))
n=0
do i=1,nfit
   if(serr(i,2).ne.0.0)then
      n=n+1
      solin(n)=sol(i) !contains subset of sol that is passed to ldif
   endif
enddo

tol=1.0d-8 !tolerance parameter for fits

!show a subset of the parameters prior to fitting

if ((sol(2).eq.0).and.(sol(3).eq.0))then
write(0,'(A130)') "   RHO        NL3        NL4        ZPT        EP1        PE1        &
BB1        RD1        EC1        ES1        KR1        TE1  "
write(0,503) sol(1),sol(4),sol(5),sol(8),sol(9),sol(10),sol(11),sol(12),sol(13),      &
   sol(14),sol(15),sol(16)
else
write(0,'(A130)') "   RHO        NL1        NL2        ZPT        EP1        PE1        &
BB1        RD1        EC1        ES1        KR1        TE1  "
write(0,503) sol(1),sol(2),sol(3),sol(8),sol(9),sol(10),sol(11),sol(12),sol(13),      &
   sol(14),sol(15),sol(16)
endif
503 format(28(1PE10.3,1X))

m=npt !number of data points
lwa=npt*nfit+5*npt*nfit !workspace size
allocate(iwa(nfit),wa(lwa),fvec(m))

!update pointers
nfit2 => nfit
sol2 => sol
serr2 => serr
nplanet2 => nplanet
aT2 => aT
aM2 => aM
aE2 => aE
aIT2 => aIT
ntt2 => ntt
tobs2 => tobs
omc2 => omc
dtype2 => dtype
nfrho2 => nfrho
rhoi2 => rhoi
rhoierr2 => rhoierr

!call fitter
call lmdif1(fcn,m,n,solin,fvec,tol,info,iwa,wa,lwa)
write(0,*) "info: ",info

n=0
do i=1,nfit
   if(serr(i,2).ne.0.0)then
      n=n+1
      sol(i)=solin(n)
   endif
enddo

!show a subset of the parameters after fitting
if ((sol(2).eq.0).and.(sol(3).eq.0))then
   write(0,503) sol(1),sol(4),sol(5),sol(8),sol(9),sol(10),sol(11),     &
    sol(12),sol(13),sol(14),sol(15),sol(16)
else
   write(0,503) sol(1),sol(2),sol(3),sol(8),sol(9),sol(10),sol(11),     &
    sol(12),sol(13),sol(14),sol(15),sol(16)
endif

return
end subroutine fittransitmodel3

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn(m,n,x,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittermod
implicit none
integer :: m,n,iflag
real(double) :: x(n),fvec(m)
!local vars
integer :: i,j,nfit
real(double) :: y,yy,drho,rhoin(9),yp,dsig,c1,c2,c3,c4
real(double), allocatable, dimension(:) :: sol
data rhoin/-1.0d2,-3.0d0,-2.0d0,-1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,1.0d2/

nfit=nfit2
allocate(sol(nfit))
sol(1:nfit)=sol2(1:nfit)

j=0
do i=1,nfit
   if(serr2(i,2).ne.0.0)then
      j=j+1
      sol(i)=x(j)
   endif
   if(j.gt.n)write(0,*) "whoops.. j>n"
enddo

!Check that limb-darkening is valid
c1=sol(2)      !limb-darkening
c2=sol(3)
c3=sol(4)
c4=sol(5)
if((c1.eq.0.0).and.(c2.eq.0.0))then
   if((c3.lt.0.0).or.(c3.gt.1.0).or.(c4.lt.0.0).or.(c4.gt.1.0))then
      fvec=9.9e30
      return
   endif
elseif((c3.eq.0.0).and.(c4.eq.0.0))then
   !write(0,*) "Quad..",c1,c2
   if((c1+c2.gt.1.0).or.(c1.lt.0).or.(c1+2.0d0*c2.lt.0))then
      !write(0,*) "invalid.."
      fvec=9.9e30
      return
   endif
endif

call transitmodel(nfit,nplanet2,nplanet2,sol,m,m,aT2,aIT2,ntt2,tobs2,   &
  omc2,fvec,dtype2)
!write(0,*) aM2(1),fvec(1)

!fold in priors on rho-star
if(nfrho2.eq.0)then
   y=0.0d0
   yy=0.0d0
   do i=1,m
      y=y+(fvec(i)-aM2(i))/aE2(i)
      yy=yy+y*y
   enddo
   drho=1.0d3*sol(1)-rhoi2
   call getrhosig(rhoierr2,rhoin,9,drho,dsig)
   yp=sqrt( (yy+dsig*dsig)/yy )
else
   yp=1.0d0
endif

do i=1,m
   fvec(i)=(fvec(i)-aM2(i))/aE2(i)*yp
!   write(6,*) i,fvec(i)
!   read(5,*)
enddo
!write(0,*) "chi: ",sum(fvec(1:m))
!read(5,*)


return
end
