subroutine pdm(nstep,p,pow,npt,time,flux,ferr,freq1,freq2,nyq,ofac,     &
 Mstar,Rstar,bper,bpow)
use precision
implicit none
integer :: npt,nstep
real(double) :: ofac,nyq,freq1,freq2,Mstar,Rstar,bper,bpow
real(double), dimension(:) :: p,pow,time,flux,ferr
!local vars
integer :: nbin,i,bin,iter
integer, allocatable, dimension(:) :: nind,pf
real :: tstart,told,tcurrent,tfinish,dt,hour,minute,second
real(double) :: steps,df0,f,q,df,perold,perdone,rbin,chisq,lconst,logl, &
 rnstep
real(double), allocatable, dimension(:) :: phase,median
real(double), allocatable, dimension(:,:) :: binnedflux,binnedferr
character(80) :: tmpc

steps=ofac*(freq2-freq1)*npt/nyq
df0=(freq2-freq1)/steps

allocate(phase(npt),pf(npt))

!constant for calculate of log-likelihood
lconst=dble(npt)*0.798179868d0+sum(log(ferr*ferr))

f=freq1

!Track and estimate CPU time.
CALL CPU_TIME(tstart)
told=tstart
perold=0.0
perdone=0.0

rnstep=dble(nstep)
bpow=-99.9e30

iter=0
do while ((f.lt.freq2).and.(iter.le.nstep)) !loop over all frequencies
   iter=iter+1  !count number of steps
   p(iter)=1.0d0/f                     !search period
   call qfunc(f,Mstar,Rstar,q) !optimum duration
   nbin=int(2.0d0/q)     !complexity of period search
   rbin=dble(nbin)

   phase=time/p(iter)-floor(time/p(iter))  !Calculate phase

   !break out data according to bin
   allocate(binnedflux(npt,nbin),nind(nbin),median(nbin),               &
    binnedferr(npt,nbin))
   nind=0 !number of data points in each bin.
   do i=1,npt
      bin=int(phase(i)*rbin)+1
      if((bin.le.nbin).and.(bin.gt.0))then
         nind(bin)=nind(bin)+1
         binnedflux(nind(bin),bin)=flux(i)
         binnedferr(nind(bin),bin)=ferr(i)
      endif
   enddo

   !calculate medians, chi-sq
   chisq=0.0d0
   do i=1,nbin
      if(nind(i).gt.3)then
         call rqsort(nind(i),binnedflux(1:nind(i),i),pf)
         median(i)=binnedflux(pf(nind(i)/2),i)
         chisq=chisq+sum( ((binnedflux(1:nind(i),i)-median(i))/         &
          binnedferr(1:nind(i),i) )**2.0)
      endif
   enddo
   deallocate(binnedflux,nind,median,binnedferr)

   logl=-0.5d0*(lconst+chisq)
   !pow(iter)=-2.0d0*logl+rbin*log(dble(npt))
   pow(iter)=logl

   if(pow(iter).gt.bpow)then
      bper=p(iter)
      bpow=pow(iter)
   endif

   !write(6,*) p(iter),pow(iter)

   CALL CPU_TIME(tcurrent)
   if(tcurrent-told.gt.1.0)then
      perdone=dble(iter)/rnstep*100.0
      dt=tcurrent-tstart
      tfinish=dt/perdone*100.0-dt
      hour=tfinish/3600.0
      minute=(tfinish-floor(hour)*3600.0)/60.0
      second=(tfinish-floor(hour)*3600.0-floor(minute)*60.0)
      perold=perdone  !update 'old' values
      told=tcurrent
      write(tmpc,'(" Percent Done: ",F6.2," ETA:",I3,"h",I2.2,"m",I2.2,"s")')&
       perold,int(hour),int(minute),int(second)
      call ovrwrt(tmpc,2)
   endif

   !Update frequency for next step.
   df=q*df0
   f=f+df


enddo
write(0,*)

return
end subroutine pdm
