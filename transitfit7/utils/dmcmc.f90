subroutine dmcmc(n,fitlookup,nbodies,npt,tol,y,yserr,sol,serr,time,    &
   flux,ferr,exptime,ans,ochi,nupcor,seed,ng,flag,nbuffer,buffer,nmov, &
   gscale,corscale,nup,nb,nupdate,npars,ngcor,ngcorsub,ngprob,         &
   ngprobsub,nacor,nacorsub,naprob,naprobsub,nas,ngs,bchi)
use precision
implicit none
type (fitlookuptable), dimension(:) :: fitlookup
integer :: nbodies,npt,nupcor,i,n,seed,ng,flag,nbuffer,nmov,j,         &
   nup,nb,nupdate,npars,nacor,nacorsub,naprob,naprobsub,nas
integer, dimension(:) :: ngcor,ngcorsub,ngprob,ngprobsub,ngs
real(double) :: tol,ochi,chi1,chi2,ran2,gasdev,rchi,fratio,u,alpha,    &
   mcmctype,nsel,nsel2,corscale,acorsub,gcorsub,bchi
real(double), dimension(:) :: y,sol,time,exptime,ans,flux,ferr,gscale
real(double), allocatable, dimension(:) :: y2,sol2,m2,y3,gratio,merr
real(double), allocatable, dimension(:,:) :: mserr
real(double), dimension(:,:) :: yserr,serr,buffer
character(80) :: cout

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

flag=0

nupcor=0 !default we do not update Gibbs samplers or Correlations

chi1=ochi !re-use previous chi-square

allocate(sol2(size(sol)),m2(size(m)),y2(size(y)))
sol2=sol !make copies for sampling
m2=m
y2=y

ng=int(ran2(seed)*dble(n)+1.0d0) !choose element to change
if((ng.lt.1).or.(ng.gt.n))then !make sure we have valid choice
   return
endif

mcmctype=ran2(seed)

if((nmov.lt.nbuffer).or.(mcmctype.le.0.5))then
   select case(fitlookup(ng)%fittype) !use lookup table to update parameter
      case(1)
         sol2(fitlookup(ng)%fitparm)=gasdev(seed)*  &
            serr(fitlookup(ng)%fitparm,2)+sol(fitlookup(ng)%fitparm)*gscale(ng)
      case(2)
         m(fitlookup(ng)%fitparm)=gasdev(seed)*     &
            mserr(fitlookup(ng)%fitparm,2)+m2(fitlookup(ng)%fitparm)*gscale(ng)
      case(3)
         y2(fitlookup(ng)%fitparm)=gasdev(seed)*    &
            yserr(fitlookup(ng)%fitparm,2)+y(fitlookup(ng)%fitparm)*gscale(ng)
   end select
else
   nsel=int(ran2(seed)*dble(nbuffer-1)+1.0d0)
   nsel2=int(ran2(seed)*dble(nbuffer-1)+1.0d0)
   do i=1,n
      select case(fitlookup(ng)%fittype) !use lookup table to update parameter
         case(1)
            sol2(fitlookup(ng)%fitparm)=(buffer(i,nsel2)-buffer(i,nsel))*corscale
         case(2)
            m(fitlookup(ng)%fitparm)=(buffer(i,nsel2)-buffer(i,nsel))*corscale
         case(3)
            y2(fitlookup(ng)%fitparm)=(buffer(i,nsel2)-buffer(i,nsel))*corscale
      end select
   enddo
endif

!y input for lcmodel gets altered
allocate(y3(size(y2)))
y3=y2 !make a copy to destroy
call lcmodel(nbodies,npt,tol,y3,sol2,time,exptime,ans)
chi2=0.0d0
do i=1,npt  !calculate chi-square
   chi2=chi2+(flux(i)-ans(i))*(flux(i)-ans(i))/(ferr(i)*ferr(i))
enddo

open(unit=10,file="lctest.dat")
do i=1,npt
   write(10,503) time(i),flux(i),ans(i)
enddo
503 format(5(1X,1PE17.10))
close(10)

fratio=exp(bchi*0.5d0*(chi1-chi2))
u=ran2(seed)
alpha=min(fratio,1.0d0)
if(u.le.alpha)then
   flag=0  !accept chain
else
   flag=1  !reject chain
endif

if(sol2(1).lt.0.0d0) flag=1 !reject chains with rhostar < 0

select case(flag)
   case(0)
      ochi=chi2 !update current chi-squared
      sol=sol2 !use new solution.  Masses in common block can be left alone
      y=y2
   case default
      ochi=chi1 !keep current chi-squared as previous solution
      m=m2 !restore old masses to common block
end select

!write(0,*) "Chi: ",chi1,chi2

!Filling buffer
if((flag.eq.0).and.(nmov.lt.nbuffer))then
   nmov=nmov+1
!   write(0,*) "nmov: ",nmov
   do i=1,n
      select case(fitlookup(ng)%fittype) !use lookup table to update parameter
         case(1)
            buffer(i,nmov)=sol2(fitlookup(ng)%fitparm)
         case(2)
            buffer(i,nmov)=m(fitlookup(ng)%fitparm)
         case(3)
            buffer(i,nmov)=y2(fitlookup(ng)%fitparm)
      end select
   enddo
   write(cout,501) "nmov:",nmov
   501 format(A5,I4)
   call ovrwrt(cout,2)
endif

!Update buffer if accepted
if((flag.eq.0).and.(nup.ge.npars).and.(nmov.eq.nbuffer).and.           &
  (mcmctype.le.0.5))then
   nup=0
   nb=nb+1
   if(nb.eq.nbuffer)then
      nupcor=1
      nupdate=nupdate+1
      npars=npars+2  !keep a larger history..
   endif
   if(nb.gt.nbuffer) nb=1
   write(cout,502) "nb: ",nb,nupdate
   502 format(A4,I4,1X,I4)
   call ovrwrt(cout,2)
   do i=1,n
      select case(fitlookup(ng)%fittype) !use lookup table to update parameter
         case(1)
            buffer(i,nb)=sol2(fitlookup(ng)%fitparm)
         case(2)
            buffer(i,nb)=m(fitlookup(ng)%fitparm)
         case(3)
            buffer(i,nb)=y2(fitlookup(ng)%fitparm)
      end select
   enddo
elseif((flag.eq.0).and.(nup.lt.npars).and.(nmov.eq.nbuffer).and.       &
  (mcmctype.le.0.5))then
   nup=nup+1
endif

!scale counts
if((nmov.lt.nbuffer).or.(mcmctype.le.0.5))then
   if(flag.eq.0)then
      ngcor(ng)=ngcor(ng)+1
      if(nupdate.gt.0) ngcorsub(ng)=ngcorsub(ng)+1
   endif
   ngprob(ng)=ngprob(ng)+1
   if(nupdate.gt.0) ngprobsub(ng)=ngprobsub(ng)+1
else
   if(flag.eq.0)then
      nacor=nacor+1
      if(nupdate.gt.0) then
         nacorsub=nacorsub+1
!         write(0,*) ng,ng2,nacorsub(ng,ng2)
      endif
   endif
   naprob=naprob+1
   if(nupdate.gt.0) then
      naprobsub=naprobsub+1
!      write(0,*) ng,ng2,naprobsub(ng,ng2)
   endif
endif

!Update corscale and gscale
if((nupcor.eq.1).and.(nupdate.gt.1))then
   write(0,*) "Updating corscale:",nupdate-1
!   if(naprob.ne.naprobsub)then
   if((naprob.ne.naprobsub).and.(nas.eq.0))then
      write(0,*) "nas: ",real(nacorsub)/real(naprobsub)
      if((real(nacorsub)/real(naprobsub).lt.0.2).or.(real(nacorsub)/   &
       real(naprobsub).gt.0.3))then
         acorsub=dble(nacor-nacorsub)/dble(naprob-naprobsub)
         corscale=corscale*(0.75*(acorsub+0.01)/(0.25*                 &
            (1.0d0-acorsub+0.01)))**0.25d0
      else
         nas=1 !stop updating scale factor
      endif
   endif
   write(0,*) "corscale"
   write(0,*) nacorsub,naprobsub,dble(nacorsub)/dble(naprobsub)
   nacorsub=0 !reset counter
   naprobsub=0 !reset counter

   allocate(gratio(n))
   do i=1,n
!      if(ngprob(i).ne.ngprobsub(i))then
      if((ngprob(i).ne.ngprobsub(i)).and.(ngs(i).eq.0))then
         write(0,*) "ngs:",i,real(ngcorsub(i))/real(ngprobsub(i))
         if((real(ngcorsub(i))/real(ngprobsub(i)).lt.0.2).or.          &
          (real(ngcorsub(i))/real(ngprobsub(i)).gt.0.3))then
            gcorsub=dble(ngcor(i)-ngcorsub(i))/                        &
             dble(ngprob(i)-ngprobsub(i))
            gscale(i)=gscale(i)*(0.75*(gcorsub+0.01)/(0.25*            &
             (1.0d0-gcorsub+0.01)))**0.25d0
         else
            ngs(i)=1 !stop updating scaling factor
         endif
      endif
!     write(0,*) 'ng:',ngcorsub(i),ngprobsub(i)
      gratio(i)=dble(ngcorsub(i))/dble(ngprobsub(i))
      ngcorsub(i)=0 !reset counter
      ngprobsub(i)=0 !reset counter
   enddo
!   write(0,500) corscale,acorsub
   500 format(20(F5.2,1X))
   write(0,*) "gscale"
   write(0,500) (gscale(i),i=1,18)
   write(0,500) (gratio(i),i=1,18)
   write(0,*)
   nupcor=0 !done updating
   deallocate(gratio) !don't need this variable anymore.
endif

return
end subroutine dmcmc
