program transitfind3
use precision
implicit none
integer :: iargc,filestatus,nunit,nmax,npt,nb,in1,in2,nstep,nplot
real(double) :: nyq,mint,maxt,fread,freq1,freq2,ofac,Mstar,Rstar,qtran, &
   depth,bper,bpow,epo,pmean,std,SNR,tdur,stat1,stat2,phc
real(double), allocatable, dimension(:) :: time,flux,ferr,psmooth,p,    &
   freqs,sn
character(80) :: obsfile,cline

interface
   subroutine readkeplc(nunit,nmax,npt,time,flux,ferr)
      use precision
      integer :: nunit
      integer :: nmax
      integer :: npt
      real(double), dimension(:) :: time,flux,ferr
   end subroutine readkeplc
   subroutine nyquest(npt,time,nyq)
      use precision
      integer :: npt
      real(double), dimension(:) :: time
      real(double) :: nyq
   end subroutine nyquest
   subroutine bls2(Mstar,Rstar,freq1,freq2,nyq,ofac,nb,npt,time,flux,   &
      in1,in2,qtran,depth,bper,bpow,freqs,p,psmooth,sn)
      use precision
      integer :: npt,nb
      integer :: in1,in2
      real(double) :: Mstar,Rstar,freq1,freq2,nyq,ofac
      real(double) :: qtran,depth,bper,bpow
      real(double), dimension(:) :: time,flux,psmooth,sn
      real(double), dimension(:) :: freqs,p
   end subroutine bls2
   subroutine calcnsteps(nstep,Mstar,Rstar,freq1,freq2,nyq,ofac,n)
      use precision
      integer :: nstep
      integer :: n
      real(double) :: Mstar,Rstar,freq1,freq2,nyq,ofac
   end subroutine calcnsteps
   subroutine smooth(nstep,p,psmooth,sn,ofac)
      use precision
      integer :: nstep
      real(double), dimension(:) :: p
      real(double), dimension(:) :: psmooth,sn
      real(double) :: ofac
   end subroutine smooth
   subroutine stats(epo,bper,tdur,npt,time,flux,pmean,std,depth,snr)
      use precision
      integer :: npt
      real(double) :: epo,bper,tdur
      real(double), dimension(:) :: time,flux
      real(double) :: pmean,std,depth,snr
   end subroutine stats
   subroutine makeplot(nplot,nstep,freqs,p,obsfile,epo,bper,bpow,depth, &
    pmean,std,snr,npt,time,flux,tdur,stat1,stat2)
      use precision
      integer :: nplot,nstep,npt
      real(double) :: bper,epo,bpow,depth,pmean,std,snr,tdur,&
              stat1,stat2
      real(double), dimension(:) :: freqs,p,time,flux
      character(80) :: obsfile
   end subroutine makeplot
   subroutine medtest(nb,npt,time,flux,epo,period,tdur,stat1,stat2)
      use precision
      implicit none
      integer :: npt,nb
      real(double) :: epo,period,tdur,stat1,stat2
      real(double), dimension(:) :: time,flux
   end subroutine medtest
end interface

if(iargc().lt.1) then !check number of command line arguments
   write(0,*) "Usage: transitfind3 <filename> [nplot] [freq1] [freq2] [Msun] [Rsun]"
   write(0,*) "  <filename> - photometry (time,flux)"
   write(0,*) "  [nplot] - 0: no plot, 1:xwindow 2:ps 3:png"
   write(0,*) "  [Per1] - low Period (days) to scan, -1 for default"
   write(0,*) "  [Per2] - high Period (days) to scan, -1 for default"
   write(0,*) "  [Msun]  - mass of host star (Msun) "
   write(0,*) "  [Rsun]  - radius of host star (Rsun) "
   stop
endif
call getarg(1,obsfile) !get filename for input data

!First argument gives us the photometry file. Open it and read in data.
nunit=10
open(unit=nunit,file=obsfile,iostat=filestatus,status='old')
if(filestatus>0)then
   write(0,*) "Cannot open ",obsfile
   stop
endif

nmax=1800000 !typical max-size needed for Kepler-LC data
allocate(time(nmax),flux(nmax),ferr(nmax)) !allocate arrays

call readkeplc(nunit,nmax,npt,time,flux,ferr) !read in photometry
close(nunit)!release unit number as we are done with the file.

if(npt.lt.3)then  !we should have at least 3 data points.
   bper=0.0
   epo=0.0
   bpow=0.0
   SNR=0.0
   qtran=0.0
   depth=0.0

else !if we have data, then do all the work..

   call nyquest(npt,time,nyq) !get nyquest frequency
   !write(0,*) "nyq: ",nyq,npt

   mint=minval(time(1:npt)) !get length of observations (days)
   maxt=maxval(time(1:npt))
   freq1=2.0/(maxt-mint) !set lowest frequency
   if(iargc().ge.4)then   !check if commandline wants a different value
      call getarg(4,cline)
      read(cline,*) fread
      if(fread.gt.0.0d0) freq1=1.0/fread
   endif

   !freq2=2.0 !0.5 day min period
   freq2=(1.0/2.0)+0.01 !3.0 day min period
   freq2=0.45
   !write(0,*) "F:",freq1,freq2
   if(iargc().ge.3)then
      call getarg(3,cline)
      read(cline,*) fread
      if(fread.gt.0.0d0) freq2=1.0/fread
   endif

   write(0,*) "freqs: ",freq1,freq2

   ofac=8.0 !over sampling

   Mstar=1.0  !Mass and radius of star - to be read in
   if(iargc().ge.5)then   !check if commandline wants a different value
      call getarg(5,cline)
      read(cline,*) fread
      if(fread.gt.0.0d0) Mstar=fread
   endif
   Rstar=1.0
   if(iargc().ge.6)then   !check if commandline wants a different value
      call getarg(6,cline)
      read(cline,*) fread
      if(fread.gt.0.0d0) Rstar=fread
   endif
   write(0,*) "Mass,Radius of host: ",Mstar,Rstar
   nb=2000
   write(0,*) "Calc number of steps"
   call calcnsteps(nstep,Mstar,Rstar,freq1,freq2,nyq,ofac,npt) !number of f's
   write(0,*) "nstep: ",nstep
   allocate(psmooth(nstep),p(nstep),freqs(nstep),sn(nstep)) !allocate arrays
   psmooth=0.0d0 !intitalize to zero (1/f estimate)
   sn=1.0d0
   write(0,*) "BLS pass 1"
   call bls2(Mstar,Rstar,freq1,freq2,nyq,ofac,nb,npt,time,flux,in1,in2, &
      qtran,depth,bper,bpow,freqs,p,psmooth,sn)
   write(0,*) "Smooth"
   call smooth(nstep,p,psmooth,sn,ofac) !used to get rid of 1/f ramp
   write(0,*) "BLS pass 2"
   call bls2(Mstar,Rstar,freq1,freq2,nyq,ofac,nb,npt,time,flux,in1,in2, &
      qtran,depth,bper,bpow,freqs,p,psmooth,sn) !run BLS again but with 1/f estimate

   !write(0,*) mint,bper,in1,in2,nb

   !pad in1,in2 by 1 and make sure we don't step outside [1,nb]
   in1=in1-1
   if(in1.lt.1) in1=in1+nb
   in2=in2+1
   if(in2.gt.nb) in2=in2-nb
   if(in1.eq.in2) in2=in1+1

   !write(0,*) mint,bper,in1,in2,nb

   !Time of first Transit
   !if(in1.lt.in2)then
   !   epo=mint+bper*dble((in1+in2)/2.0)/dble(nb)
   !else
   !   epo=mint+bper*dble((in1+in2-nb)/2.0)/dble(nb)
   !endif

   if(in1.lt.in2)then
           phc = dble((in1+in2)/2.0)/dble(nb)
   else
           phc = dble((in1+in2-nb)/2.0)/dble(nb)
   endif
   epo = bper * dble( int(mint/bper)+1 ) + phc * bper

   !Transit duration.  (units = period = days)
   if(in1.le.in2)then
      tdur=bper*dble(in2-in1)/dble(nb)
   else
      tdur=bper*dble(in2+nb-in1)/dble(nb)
   endif

   !get some stats about the potential event.
   !write(0,*) "Calculating Stats"
   call stats(epo,bper,tdur,npt,time,flux,pmean,std,depth,snr)

endif

call medtest(nb,npt,time,flux,epo,bper,tdur,stat1,stat2)

!out period,T0,bpow,S/N,Tdur
bper=max(1.0e-10,bper)
write(6,500) bper,epo,bpow,SNR,tdur,depth,stat1,stat2
500  format(8(1PE17.10,1X))

nplot=0 !disable plotting by default
if(iargc().ge.2)then
   call getarg(2,cline)
   read(cline,*) nplot
endif
if(nplot.ge.1)then
   call makeplot(nplot,nstep,freqs,p,obsfile,epo,bper,bpow,depth,pmean, &
      std,snr,npt,time,flux,tdur,stat1,stat2)
endif

end program transitfind3

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine medtest(nb,npt,time,flux,epo,period,tdur,stat1,stat2)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: npt,nb
real(double) :: epo,period,tdur,stat1,stat2
real(double), dimension(:) :: time,flux
!local vars
integer :: i,j,k,midx
integer, dimension(:), allocatable :: p
real(double) :: t,dnb,poff,phase,ph0,phs,phe,m,maxm,event0,event1,event2, &
   ph1,ph2,dph,stdev,pmean,std
real(double), dimension(:), allocatable :: f,med,ph

allocate(f(max(nb,npt)),p(max(nb,npt)),med(nb),ph(nb))


phs=-tdur/period/2.0d0 !- 1/2 tdur - these are relative to zero.                                 
phe=+tdur/period/2.0d0 !+ 1/2 tdur

dnb=dble(nb)
do i=1,nb
   ph0=dble(i-1)/dnb
   ph(i)=ph0

   poff=ph0*period
   k=0
   do j=1,npt
      !phase data
      t=time(j)-(epo+poff)
      phase=t/period-int(t/period) !Transit is located at phase=0
      if (phase.lt.0.0) phase=phase+1.0d0 !correct for negative phase
      if (phase.gt.0.5) phase=phase-1.0d0 !phase goes from -0.5 to 0.5
      if ((phase.gt.phs).and.(phase.lt.phe))then
         k=k+1
         f(k)=flux(j)
      endif
   enddo

   if(k.gt.2)then
      call rqsort(k,f,p)
      med(i)=-f(p((k)/2))
   else
      med(i)=0.0d0
   endif
   !write(0,*) ph0,med(i)

enddo

event0=med(1) !this is the primary event

!get next event
ph0=0.0
phs=ph0-1.5*tdur/period+1.0 !-tdur - these are relative to zero.
phe=ph0+1.5*tdur/period !+tdur
maxm=0.0
midx=1
do i=1,nb
   if((ph(i).gt.phe).and.(ph(i).lt.phs))then
      m=abs(med(i))
      if(m.gt.maxm)then
        midx=i
        maxm=m
      endif
   endif
enddo

!write(0,*) "maxm: ",midx,med(midx)
event1=med(midx)
ph1=ph(midx)

!find third event
dph=1.5*tdur/period
maxm=0.0
midx=1
k=0
do i=1,nb
   if((ph(i).gt.phe).and.(ph(i).lt.phs).and.(abs(ph(i)-ph1).gt.dph))then
      m=abs(med(i))
      if(m.gt.maxm)then
        midx=i
        maxm=m
      endif
      k=k+1
      f(k)=med(i)
   endif
enddo
std=stdev(k,f,pmean) !get Standard deviation for stats

!write(0,*) "maxm: ",midx,med(midx)
event2=med(midx)
ph2=ph(midx)

!write(0,*) event0,event1,event2,std

stat1=(abs(event0)-abs(event1))/std
stat2=(abs(event0)-abs(event2))/std

!write(0,*) "Diff-med test: ",stat1,stat2

deallocate(f,p)

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine makeplot(nplot,nstep,freqs,p,obsfile,epo,bper,bpow,depth,    &
   pmean,std,snr,npt,time,flux,tdur,stat1,stat2)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nplot,nstep,npt
real(double) :: epo,bper,phase,bpow,depth,pmean,std,snr,tdur,stat1,stat2
real(double), dimension(:) :: freqs,p,time,flux
real(double), allocatable, dimension(:) :: pers
character(80) :: obsfile,txtout

select case(nplot)
   case(1)
      call pgopen('/xserve')
   case(2)
      call pgopen('find.ps/vcps')
   case(3)
      call pgopen('find.png/png')
end select

call pgsch(2.9)
call pgsubp(1,4)
call pgvport(0.1,0.95,0.2,0.8) !gives enough room for labels

call pgpage()

if(npt.ge.3)then

   allocate(pers(nstep))
   pers=1.0d0/freqs
   call plot(nstep,pers,p,obsfile,"Period (d)","P",0,0.0d0,0.0d0)
   call pgpage()
   call plot(npt,time,flux,obsfile,"BJD-2454900","Flux",1,epo,bper)
   call pgpage()
   phase=epo/bper-int(epo/bper)
   if(phase.lt.0.0d0) phase=phase+1.0d0
   call plotph(npt,time,flux,bper,phase)
   call pgpage()
   call pgwindow(0.0,1.0,0.0,1.0)

   write(txtout,501) "Per: ",bper
   501 format(A5,1X,1PE17.10)
   call pgptxt(0.1,0.9,0.0,0.5,txtout)
   write(txtout,501) "Epo: ",epo
   call pgptxt(0.1,0.77,0.0,0.5,txtout)
   write(txtout,501) "Pow: ",bpow
   call pgptxt(0.1,0.64,0.0,0.5,txtout)
   write(txtout,501) "S/N: ",snr
   call pgptxt(0.1,0.51,0.0,0.5,txtout)
   write(txtout,501) "Dur: ",tdur
   call pgptxt(0.1,0.38,0.0,0.5,txtout)
   write(txtout,501) "Dep: ",depth*1.0e6
   call pgptxt(0.1,0.25,0.0,0.5,txtout)
   write(txtout,501) "d12: ",stat1
   call pgptxt(0.1,0.12,0.0,0.5,txtout)
   write(txtout,501) "d13: ",stat2
   call pgptxt(0.1,-0.01,0.0,0.5,txtout)

   call plottrans2(npt,time,flux,bper,phase,tdur)

endif

call pgclos()

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine stats(epo,period,tdur,npt,time,flux,fmean,std,depth,snr)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: npt
real(double) :: epo,period,tdur,fmean,std,depth,snr
real(double), dimension(:) :: time,flux
!local vars
integer :: i,k,kin
integer, dimension(:), allocatable :: p
real(double) :: tdurfac,ph1,ph2,toff,phase,stdev
real(double), dimension(:), allocatable :: f,fin

allocate(f(npt),fin(npt),p(npt))

tdurfac=0.5d0 !exactly 1 transit duration

ph1=0.75-tdurfac*tdur/period
if(ph1.lt.0.5)ph1=0.5
ph2=0.75+tdurfac*tdur/period
if(ph2.gt.1.0)ph2=1.0

toff=0.75-(epo/period-int(epo/period))
if(toff.lt.0.0)toff=toff+1.0

!write(0,*) 'ph1,ph2: ',ph1,ph2

k=0 !number of data points of out transit
kin=0 !number of data point in transit
do i=1,npt
   phase=time(i)/period-int(time(i)/period)
   phase=phase+toff
   if(phase.lt.0.0) phase=phase+1.0
   if(phase.gt.1.0) phase=phase-1.0
   if((phase.lt.ph1).or.(phase.gt.ph2))then
      k=k+1
      f(k)=flux(i)
   else !data is in-transit
      kin=kin+1
      fin(kin)=flux(i)
   endif
enddo

!write(0,*) 'k,kin', k,kin

!get the standard deviation of out-of-transit data
if(k.gt.3)then
   std=stdev(k, f, fmean) !fmean is the mean of out-out-transit data.  
else
   fmean=0.0
   std=0.0
endif

if(kin.gt.1)then
   !get the median of in-transit data
   call rqsort(kin,fin,p)
   depth=-fin(p((kin)/2))+fmean !by definition, transits have a positive depth
else
   depth=fmean
endif

if((std.gt.0.0).and.(kin.gt.0))then
   snr=depth/std*sqrt(dble(kin))
else
   snr=0.0
endif

!write(0,*) "stats: ",tdur,std,fmean,meddepth

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine smooth(nstep,p,psmooth,sn,ofac)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nstep,i,i1,i2,k,j,nsamp
integer, allocatable, dimension(:) :: nd
real(double) :: ofac,pmean,stdev
real(double), dimension(:) :: p,psmooth,sn
real(double), allocatable, dimension(:) :: pcut

!nsamp=max(int(ofac),nstep/1000)
!nsamp=min(int(ofac)*100,nstep) ! original smoothing call
nsamp=min(int(ofac)*100,nstep)

allocate(pcut(nstep),nd(nstep))

do i=1,nstep
   i1=max(1,i-nsamp)
   i2=min(nstep,i+nsamp)
   k=0
   do j=i1,i2
      if((j.gt.i+nsamp/4).or.(j.lt.i-nsamp/4))then !chop out frequency of interest
         k=k+1
         pcut(k)=p(j)
      endif
   enddo
!   psmooth(i)=pcut(1)
!   do j=2,k
!      psmooth(i)=min(psmooth(i),pcut(j))
!   enddo
  call rqsort(k,pcut,nd)
  psmooth(i)=pcut(nd(k/2))
  sn(i)=stdev(k, pcut, pmean)
enddo

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine calcnsteps(nstep,Mstar,Rstar,freq1,freq2,nyq,ofac,n)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer nstep,n
real(double) :: Mstar,Rstar,freq1,freq2,nyq,ofac,df0,f,q,df,steps

interface
   subroutine qfunc(f,Mstar,Rstar,q)
      use precision
      real(double) :: f,Mstar,Rstar,q
   end subroutine qfunc
end interface

!write(0,*) ofac,n,nyq
steps=ofac*(freq2-freq1)*n/nyq
df0=(freq2-freq1)/dble(steps)
f=freq1
nstep=0
!write(0,*) "steps: ",steps,freq1,freq2
do while(f.lt.freq2)
!   write(0,*) "nstep: ",nstep,f,df0
   nstep=nstep+1 !increase counter to count number of steps
   call qfunc(f,Mstar,Rstar,q)
   df=q*df0 !estimate next frequency step based on optimum q
   f=f+df !next frequency to scan
enddo

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine bls2(Mstar,Rstar,freq1,freq2,nyq,ofac,nb,n,t,x,in1,in2,qtran,&
   depth,bper,bpow,freqs,p,psmooth,sn)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: n,nstep,nb,kmi,kma,kkmi,minbin,nb1,nbkma,i,j,jnb,k,kk,nb2,   &
   jn1,jn2,in1,in2
integer, allocatable, dimension(:) :: ibi
real(double) :: Mstar,Rstar,freq1,freq2,f,df,q,nyq,ofac,steps,df0,rn,   &
   rnb,bpow,onehour,qma,qmi,s,t1,f0,p0,ph,power,rn1,pow,rn3,s3,sn1,ps1, &
   qtran,depth,bper
real(double), dimension(:) :: t,x,freqs,p,psmooth,sn
real(double), allocatable, dimension(:) :: y,u,v

interface
   subroutine qfunc(f,Mstar,Rstar,q)
      use precision
      real(double) :: f,Mstar,Rstar,q
   end subroutine qfunc
end interface

!The BLS code is buggy and ibi will segfault with overflows if you just
!allocate nb worth of space.  So I give twice as much space as asked for
allocate(y(nb*2),ibi(nb*2),u(size(t)),v(size(t)))

steps=ofac*(freq2-freq1)*n/nyq
df0=(freq2-freq1)/dble(steps)
!write(0,*) "x:",n,steps,nyq

!write(0,*) "x:",freq1,freq2,df0
f=freq1
nstep=0
rn=dfloat(n)
rnb=dfloat(nb)
bpow=0.0d0
onehour=1.0/24.0 !one hour convected to day.
minbin=5
nb1   = nb+1

!
!=================================
!     Set temporal time series
!=================================
!
s=0.0d0
!t1=t(1)
t1=minval(t)
do i=1,n
   u(i)=t(i)-t1
   s=s+x(i)
enddo
s=s/rn
do i=1,n
   v(i)=x(i)-s
enddo

do while(f.lt.freq2)
   nstep=nstep+1 !increase counter to count number of steps
   f0=f  !BLS cuts uses f0, so just copy
   p0=1.0d0/f0  !period
   freqs(nstep)=f0

!  Added by JR
!   sn1=sn(nstep)
   sn1=1.0 !K
   ps1=psmooth(nstep)
!   sn1=psmooth(nstep)
!   ps1=0.0!
!  Added by JR

!  optimum q to scan based on stellar parameters
   call qfunc(f,Mstar,Rstar,q)
!  range of q to scan (right now it's q/2 and 2*q and no shorter than 1 hour
   qmi=max(onehour*f,q/2.0) !don't search for durations less than 1 hour
   qma=min(0.05,q*2.0) !don't search for durations greater that 5% of period
   if(qma.lt.qmi) qma=min(qmi*2.0,0.4) !make sure qma > qmi

!  BLS guts
   kmi=idint(qmi*rnb)
   if(kmi.lt.1) kmi=1
   kma=idint(qma*rnb)+1
   kkmi=idint(rn*qmi)
   if(kkmi.lt.minbin) kkmi=minbin
   nbkma = nb+kma

!
!======================================================
!     Compute folded time series with  *p0*  period
!======================================================
!
   do j=1,nb
      y(j) = 0.0d0
      ibi(j) = 0
   enddo

   do i=1,n
      ph     = u(i)*f0
      ph     = ph-idint(ph)
      j      = 1 + idint(nb*ph)
      ibi(j) = ibi(j) + 1
      y(j) =   y(j) + v(i)
   enddo

!
!-----------------------------------------------
!     Extend the arrays  ibi()  and  y() beyond
!     nb   by  wrapping
!
   do j=nb1,nbkma
      jnb    = j-nb
      ibi(j) = ibi(jnb)
        y(j) =   y(jnb)
   enddo
!-----------------------------------------------
!
!===============================================
!     Compute BLS statistics for this period
!===============================================
!
   power=0.0d0
!
   do i=1,nb
      s     = 0.0d0
      k     = 0
      kk    = 0
      nb2   = i+kma
      do j=i,nb2
         k     = k+1
         kk    = kk+ibi(j)
         s     = s+y(j)
         if(k.lt.kmi) cycle
         if(kk.lt.kkmi) cycle
         rn1   = dfloat(kk)
         pow   = s*s/(rn1*(rn-rn1))
         if(pow.lt.power) cycle
         power = pow
         jn1   = i
         jn2   = j
         rn3   = rn1
         s3    = s
      enddo
   enddo

   power = dsqrt(power)
!   if(nstep.le.30) then
!      power=0.0d0
!      write(0,*) power,qmi,qma
!   endif
   p(nstep) = (power-ps1)/sn1

!   write(6,*) f0,p(nstep)

   if((power-ps1)/sn1.gt.bpow)then
      bpow  =  (power-ps1)/sn1
      in1   =  jn1
      in2   =  jn2
      qtran =  rn3/rn
      depth = -s3*rn/(rn3*(rn-rn3))
      bper  =  p0
   endif

   df=q*df0 !estimate next frequency step based on optimum q
   f=f+df !next frequency to scan
enddo

write(0,*) "steps: ",nstep
if(in2.gt.nb) in2 = in2-nb

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine qfunc(f,Mstar,Rstar,q)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
real(double) :: f,q,Mstar,Rstar,fac,G,Rsun,Msun,fsec

G=6.674d-11 !N m^2 kg^-2  Gravitation constant
fac=1.083852140278d0 !(2pi)^(2/3)/pi
Rsun=696265.0d0*1000.0d0 !m  radius of Sun
Msun=1.9891d30 !kg  mass of Sun

fsec=f/86400.0d0 !convert from c/d to c/sec

q=fac*Rstar*Rsun/(G*Mstar*Msun)**(1.0/3.0)*fsec**(2.0/3.0)

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine nyquest(npt,time,nyq)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: npt,ndt,i
integer, allocatable, dimension(:) :: p
real(double) :: nyq,mode
real(double), dimension(:) :: time
real(double), allocatable, dimension(:) :: dt

allocate(dt(npt),p(npt)) !allocate arrays

ndt=0 !initalize counter
do i=2,npt
   ndt=ndt+1
   dt(ndt)=time(i)-time(i-1) !difference between times
enddo

call rqsort(ndt,dt,p) !sort data
mode=dt(p(ndt/2)) !get median value

nyq=1.0/(2.0*mode) !calculate nyquest freq

!write(0,*) "nyq:",npt,nyq

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readkeplc(nunit,nmax,npt,time,flux,ferr)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer nunit,npt,nmax,i,filestatus
real(double), dimension(:) :: time,flux,ferr

i=1
do
   if(i.gt.nmax)then !ran out of array space
      write(0,*) "not enough array space Increase nmax"
      stop
   endif
   read(nunit,*,iostat=filestatus) time(i),flux(i),ferr(i)
   if(filestatus == 0) then
      time(i)=time(i)-54900.0d0+0.5d0 !convert to BJD-2454900
!      if(records(i,29).gt.10.0)
!      write(0,*) time(i),flux(i),ferr(i)
      i=i+1
      cycle
   elseif(filestatus == -1) then  !EOF
      exit
   else
      write(0,*) "previous time read:",time(i-1)
      write(0,*) "File Error! line:",i
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
npt=i-1
!write(0,*) "Records read: ",npt

return
end
