program transitmcmc7
use precision
implicit none
type (fitlookuptable), allocatable, dimension(:) :: fitlookup
integer iargc,npt,nbodies,nmax,i,j,k,seed,niter,n,nas,npars,nacor,     &
   nacorsub,naprob,naprobsub,nmov,nup,nb,nbupdate,naccept,nupdate,     &
   flag,ii,ng,nupcor,nbuffer,nplanet
integer, dimension(3) :: now
integer, allocatable, dimension(:) :: ngs,ngcor,ngcorsub,ngprob,       &
   ngprobsub
real(double) :: tol,dumr,ran2,gasdev,bchi,fchi,ochi,corscale,          &
   accrate
real(double), allocatable, dimension(:) :: time,flux,ferr,exptime,y,   &
   yerr,sol,err,ans,gscale,y2,merr
real(double), allocatable, dimension(:,:) :: yserr,serr,buffer,mserr
character(80) :: photfile,inputsol,command

interface
   subroutine readkeplc(photfile,npt,time,flux,ferr,exptime)
      use precision
      implicit none
      character(80), intent(in) :: photfile
      integer, intent(out) :: npt
      real(double), dimension(:), intent(inout) :: time,flux,ferr,exptime
   end subroutine readkeplc
end interface

interface
   subroutine calcnbodies(inputsol,nbodies)
      use precision
      implicit none
      character(80), intent(in) :: inputsol
      integer, intent(out) :: nbodies
   end subroutine calcnbodies
end interface

interface
   subroutine readinputsol(inputsol,nbodies,y,yserr,yerr,sol,serr,err)
      use precision
      implicit none
      character(80), intent(in) :: inputsol
      integer, intent(in) :: nbodies
      real(double), dimension(:), intent(out) :: y,yerr,sol,err
      real(double), dimension(:,:), intent(out) :: yserr,serr
   end subroutine readinputsol
end interface

interface
   subroutine lcmodel(nbodies,npt,tol,y,sol,time,itime,ans)
      use precision
      implicit none
      integer, intent(in) :: nbodies,npt
      real(double), intent(in) :: tol
      real(double), dimension(:), intent(inout) :: y,ans
      real(double), dimension(:), intent(in) :: sol,time,itime
   end subroutine lcmodel
end interface

interface
   subroutine dmcmc(n,fitlookup,nbodies,npt,tol,y,yserr,sol,serr,time, &
     flux,ferr,exptime,ans,ochi,nupcor,seed,ng,flag,nbuffer,buffer,    &
     nmov,gscale,corscale,nup,nb,nupdate,npars,ngcor,ngcorsub,ngprob,  &
     ngprobsub,nacor,nacorsub,naprob,naprobsub,nas,ngs,bchi)
      use precision
      implicit none
      type (fitlookuptable), dimension(:), intent(in) :: fitlookup
      integer, intent(in) :: nbodies,npt,n,ng,nbuffer
      integer, intent(inout) :: nupcor,seed,flag,nmov,nup,nb,nupdate,  &
         npars,nacor,nacorsub,naprob,naprobsub,nas
      integer, dimension(:), intent(inout) :: ngcor,ngcorsub,ngprob,   &
         ngprobsub,ngs
      real(double), intent(in) :: tol,bchi
      real(double), intent(inout) :: ochi,corscale
      real(double), dimension(:), intent(in) :: time,exptime,          &
         flux,ferr
      real(double), dimension(:), intent(inout) :: ans,y,sol,gscale
      real(double), dimension(:,:), intent(in) :: yserr,serr
      real(double), dimension(:,:), intent(inout) :: buffer
   end subroutine dmcmc
end interface

interface
   subroutine exportfit(nplanet,sol,serr,err,y,yerr,yserr)
      use precision, only: double
      implicit none
      integer, intent(in) :: nplanet
      real(double), dimension(:), intent(in) :: sol,err,y,yerr
      real(double), dimension(:,:), intent(in) :: serr,yserr
   end subroutine exportfit
end interface


tol=1.0d-12  !default tolerance level

if(iargc().lt.3) then
   write(0,*) "Usage: transitfit7 <photometry> <inputsol> <niter>"
   stop
endif

call getarg(1,photfile)
nmax=100000
allocate(time(nmax),flux(nmax),ferr(nmax),exptime(nmax))
call readkeplc(photfile,npt,time,flux,ferr,exptime)
write(0,*) "Number of observations read: ",npt

call getarg(2,inputsol)
call calcnbodies(inputsol,nbodies)
allocate(m(nbodies),mserr(nbodies,2),merr(nbodies),y(nbodies*6),       &
   yserr(nbodies*6,2),yerr(nbodies*6),sol(8+nbodies),                  &
   serr(8+nbodies,2),err(8+nbodies),y2(nbodies*6))
write(0,*) "nbodies: ",nbodies
call readinputsol(inputsol,nbodies,y,yserr,yerr,sol,serr,err)
do i=1,nbodies
   write(0,500) m(i),(y(j),j=6*i-5,6*i)
enddo
500  format(28(1PE10.3,1X))
write(0,*) " "

call getarg(3,command)
read(command,*) niter
if(niter.lt.1)then
   write(0,*) "niter must be greater that 0"
   stop
endif

!get chi-sq
allocate(ans(npt))
y2=y !y get altered by lcmodel
call lcmodel(nbodies,npt,tol,y2,sol,time,exptime,ans)
bchi=0.0d0
do i=1,npt
   bchi=bchi+(flux(i)-ans(i))*(flux(i)-ans(i))/(ferr(i)*ferr(i))
enddo
fchi=bchi !keep tabs on the best chi-squared
write(0,*) "bchi",bchi,npt-1
bchi=dble(npt-1)/bchi
write(0,*) "bchi_f:",bchi

open(unit=10,file="lctest.dat")
do i=1,npt
   write(10,503) time(i),flux(i),ans(i)
enddo
503 format(5(1X,1PE17.10))
close(10)
write(0,*) "wait1..."
read(5,*)

!initalize random number
call itime(now)
seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
write(0,*) "Seed: ",seed
dumr=ran2(-seed)

allocate(fitlookup(8+nbodies+nbodies+nbodies*6)) !allocate maximum size needed
!how many variables are we fitting?
n=0
do i=1,8+nbodies
   if(serr(i,2).ne.0.0d0)then
      n=n+1
      fitlookup(n)%fittype=1
      fitlookup(n)%fitparm=i
   endif
enddo
do i=1,nbodies
   if(mserr(i,2).ne.0.0d0)then
      n=n+1
      fitlookup(n)%fittype=2
      fitlookup(n)%fitparm=i
   endif
enddo
do i=1,nbodies*6
   if(yserr(i,2).ne.0.0d0)then
      n=n+1
      fitlookup(n)%fittype=3
      fitlookup(n)%fitparm=i
   endif
enddo
write(0,*) "n:",n

!allocate flag for updating jump scale
allocate(ngs(n),gscale(n),ngcor(n),ngcorsub(n),ngprob(n),ngprobsub(n))

nas=0 !initalized flag for updating jump scale
npars=1 !counting number parameters that are fitted
corscale=1.0d0
nacor=0
nacorsub=0
naprob=0
naprobsub=0

nmov=0 !initialize buffer size for dmcmc
nup=0 !when to update buffer (0=no,1=yes)
nb=0  !which buffer element to replace (oldest)
nupdate=0 !number of times scales have been updated

!initilization of counters
do i=1,n
   ngs(i)=0 ! initialize flag for updating jump scale
   gscale(i)=1.0d0
   ngcor(i)=0
   ngcorsub(i)=0
   ngprob(i)=0
   ngprobsub(i)=0
enddo

!now we giggle the best-fit values to start the MCMC.
do k=1,8+nbodies
   if(serr(k,2).ne.0.0d0)then
      sol(k)=gasdev(seed)*serr(k,2)+sol(k)
   endif
enddo
do k=1,nbodies
   if(mserr(k,2).ne.0.0d0)then
      m(k)=gasdev(seed)*mserr(k,2)+m(k)
   endif
enddo
do k=1,nbodies*6
   if(yserr(k,2).ne.0.0d0)then
      y(k)=gasdev(seed)*yserr(k,2)+y(k)
   endif
enddo

!calculate chi-square of jittered model
y2=y !y gets altered by lcmodel
call lcmodel(nbodies,npt,tol,y2,sol,time,exptime,ans)
ochi=0.0d0
do i=1,npt
   ochi=ochi+(flux(i)-ans(i))*(flux(i)-ans(i))/(ferr(i)*ferr(i))
enddo
write(0,*) "ochi: ",ochi

open(unit=10,file="lctest.dat")
do i=1,npt
   write(10,503) time(i),flux(i),ans(i)
enddo
close(10)
write(0,*) "wait2..."
read(5,*)

nbuffer=500 !buffer for storing previous MCMC chains
allocate(buffer(n,nbuffer))

naccept=0
do ii=1,niter
   call dmcmc(n,fitlookup,nbodies,npt,tol,y,yserr,sol,serr,time,flux,  &
      ferr,exptime,ans,ochi,nupcor,seed,ng,flag,nbuffer,buffer,nmov,   &
      gscale,corscale,nup,nb,nupdate,npars,ngcor,ngcorsub,ngprob,      &
      ngprobsub,nacor,nacorsub,naprob,naprobsub,nas,ngs,bchi)
   write(6,501) ochi,flag,fitlookup(ng)%fittype,                       &
      fitlookup(ng)%fitparm,                                           &
      (sol(i),i=1,8+nbodies),(m(i),i=1,nbodies),(y(i),i=1,nbodies*6)
   501 format(1PE17.10,3(1X,I2),108(1X,1PE17.10))
   if(flag.eq.0) naccept=naccept+1 !track acceptance rate

   if(ochi.lt.fchi)then  !check for new chi-squared min
      fchi=ochi
      bchi=dble(npt-1)/fchi
      write(0,*) "fchi",fchi,bchi
      write(0,*) "Exporting fit"
      nplanet=nbodies-1
      call exportfit(nplanet,sol,serr,err,y,yerr,yserr)
   endif

enddo
accrate=dble(naccept)/dble(niter)
write(0,*) "Accept Rate: ",accrate

end program transitmcmc7
