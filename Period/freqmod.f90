program freqmod
!Jason Rowe 2014/08/07 Version 1.0
use precision
implicit none
integer :: fixfreq,steps,nfit,nunit,npt,filestatus,nfitl,nstarmax,      &
   niter,i,nstar,plot,kID,nfmax,ma,nb,nc,itmax,npt2,slsteps,j,ndphase,  &
   kk,cc,jj,k,n
integer, parameter :: nmax=1650000
integer, allocatable, dimension(:) :: ndt,nd
real(double) :: keplertime,tbin,snlimit,freq1,freq2,nyq,nyquest,maxx,   &
   minx,ofac,sig,rmavg,stdev,mean,bper,btheta,sn,xcoo,ycoo,             &
   avg,dnu,sfreq1,sfreq2,nu,phase,dphase,rndphase,mfreq,mtheta,per,     &
   twopi,pi,ztime,atheta,afreq,mstdev,mstd,abin,pid2
real(double), dimension(nmax) :: time,mag,merr,exptime,time2,mag2,merr2
real(double), allocatable, dimension(:) :: dt,std,res,a,aerr,w,pers,    &
   sns,pers2,a2,aerr2,dfreqs,dthetas,a3,aerr3,dfreqserr,dthetaserr
character :: ans
character(80) :: filename,title,freqfile
common /data/ time,mag,merr,npt,nfit

fixfreq=0 !0=fit frequencies, 1=fixed frequencies
keplertime=2451545.00000000 !HJD of Kepler start time
tbin=0.0  !initialize binning of data parameters (0==do not bin)
snlimit=20.0 !default limit for determination of significant amps
freq1=-1.0 !set to 1/(maxT-minT)
freq2=-1.0 !set to nyquest
steps=0 !initialize number of steps (default is ~1200 minimum)
nfit=1  !for harmonic fits
nstarmax=300 !max num of freqs to fit.
nfmax=2000 !max number of fitted variables.
ofac=8.0 !over sampling

pi=3.141592654
twopi=2.0d0*pi
pid2=Pi/2.0d0

if(iargc().lt.1) then
   write(0,*) "Usage: freqmod keplerlc.dat"
   stop
else
   call getarg(1,filename)
endif

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then
   write(0,*) "Cannot open ",filename
   stop
endif

!read in Kepler photometry
call readkeplc(nunit,nmax,npt,time,mag,merr,exptime,Keplertime)
close(nunit)
title=" " !if you want a title on the top

!These lines convert to parts-per-thousand
mag=mag*1000.0
merr=merr*1000.0

write(0,*) "Points Read: ",npt

!get max/min times
maxx=maxval(time)
minx=minval(time)

call pgopen('/null')
!call pgask(.false.)
call pgsch(2.9)  !bigger text
call pgsubp(1,4) !4 panels
call pgvport(0.1,0.95,0.2,0.8) !gives enough room for labels
call pgpage()

allocate(dt(nmax),ndt(nmax))
nyq=nyquest(nmax,npt,time,dt,ndt)

nb=5
nc=2
if(freq1.lt.0.0) freq1=1.0/(maxx-minx)
if(freq2.lt.0.0) freq2=nyq*3.0 !super-nyquest searching
if(steps.le.0) steps=int(ofac*(freq2-freq1)*npt/nyq)
if(steps.gt.1e6)then
!   write(0,*) "reducing steps to 1e6"
   steps=1e6
endif
nfitl=3
niter=nstarmax-1
sig=3.0 !set level of sigma clipping

! sigma clipping (not using if counter is set to zero below)
do i=1,0
   call sigclip(npt,time,mag,merr,sig,rmavg)
enddo

call avgrm(npt,mag,rmavg) !remove zero-point from data.
if(tbin.gt.0.0) call bindt(npt,time,mag,merr,tbin) !not tested!!
!call boxcar(npt,time,mag,merr,boxbin)
!tbin=1.0*500.0
tbin=0.25*24.0*60.0
!call spdetrend(npt,time,mag,merr,tbin,sky) !spline detrending
!call sigclip(npt,time,mag,merr,sig,rmavg)
!call plotdata(npt,time,mag,merr,-1.0,-1.0,1,1,MOSTtime)
!write(6,*) id,xcoo,ycoo,npt
call avgrm(npt,mag,rmavg) !make sure zero point is up-to-date

allocate(std(2))
std(1)=1.0d3*stdev(npt,mag,mean)

allocate(pers(nstarmax),sns(nstarmax),res(nmax),a(nfmax),aerr(nfmax),   &
   w(nstarmax))
nstar=0
kid=0

write(0,*) freq1,freq2,steps

if(iargc().ge.2)then
   call getarg(2,freqfile)
   nunit=10
   open(unit=nunit,file=freqfile,iostat=filestatus,status='old')
   if(filestatus>0)then
      write(0,*) "Cannot open ",filename
      stop
   endif
   write(0,*) "Reading in Frequencies"
   call readfreqs(nunit,ma,nstar,a,sns,ztime)
   write(0,*) ma,nstar
   j=(ma-nstar)/nstar
   do k=1,nstar
      cc=(j+1)*(k-1)+2
      pers(k)=twopi/a(cc-1)
      w(k)=1.0/pers(k)
!      write(0,*) k,pers(k),w(k)
   enddo
   write(0,*) "Done reading frequencies"
   call shaperm(npt,nfit,time,mag,merr,res,pers,avg,a,1,3,plot,w,       &
       nstar,Keplertime)
else
   xcoo=-1.0d0
   call plotdata(npt,time,mag,merr,-1.0d0,-1.0d0,1,1,Keplertime,title)
   call jmfourw(npt,time,mag,merr,freq1,freq2,steps,bper,btheta,sn,1,2, &
    snlimit,1)

   write(0,*) "bper:",bper
   plot=0
   bper=1.0/bper !cosinefit works with Periods not frequencies
   nstar=nstar+1 !number of individual frequencies for fitting.
   pers(nstar)=bper !store best *period*
   sns(nstar)=sn !store SN of best period
   call cosinefit(pers,res,kid,xcoo,ycoo,avg,a,btheta,rmavg,plot,ma,aerr,&
     1,3,w,nstar,fixfreq)
   call windowfn(npt,nfit,time,mag,merr,bper,freq1,freq2,avg,a,nb,nc,   &
     steps,1,3,nstar,pers)
   plot=0
   call shaperm(npt,nfit,time,mag,merr,res,pers,avg,a,1,3,plot,w,nstar, &
     Keplertime)
   call jmfourw(npt,time,res,merr,freq1,freq2,steps,bper,btheta,sn,1,4, &
     snlimit,1)

   !plot=1
   !write(0,*) "Next peak?          "
   !read(5,501) ans
   !501  format(A1)
   !if(ans.eq."n")goto 24
   ans="y"  !fit everything above threshold.

   itmax=nstarmax-1 !up to nstarmax frequencies

   i=0 !counter

   if((sn.lt.snlimit).or.(i.gt.itmax)) then
      ans="n"
   else
      ans="y"
   endif

   do while(ans.eq."y")
      i=i+1
      write(0,*) "Iteration #: ",i
      write(0,*) "bfreq:",bper,sn
!      call pgpanl(1,4)
!      call pgpage()
      bper=1.0/bper
      nstar=nstar+1
      pers(nstar)=bper
      sns(nstar)=sn
!      call jmfour(npt,time,mag,per1,per2,steps,bper,btheta,1,1)
      call pgpanl(1,2)
      call pgeras()
!      call plotph(npt,time,res,bper,1,2)
      plot=1
      call cosinefit(pers,res,kID,xcoo,ycoo,avg,a,btheta,rmavg,plot,    &
       ma,aerr,1,2,w,nstar,fixfreq)
      write(0,503) i+1,"bfreq:",1.0/bper,sn
      503 format(i3,1x,A6,1X,F8.5,1X,F6.3)
!      call plotdata(npt,time,mag,merr,-1.0,-1.0,1,2,MOSTtime)
!      call plotfit(npt,nfit,time,mag,pers,avg,a,w,nstar)
      call pgpanl(1,3)
      call pgeras()
      call shaperm(npt,nfit,time,mag,merr,res,pers,avg,a,1,3,plot,w,       &
       nstar,Keplertime)
      call pgpanl(1,4)
!      call pgeras()
      call jmfourw(npt,time,res,merr,freq1,freq2,steps,bper,btheta,        &
       sn,1,4,snlimit,1)
!      call pdm(npt,time,res,nb,nc,per1,per2,steps,bper,btheta,1,4)
      if((sn.lt.snlimit).or.(i.gt.itmax)) then
         ans="n"
      else
         ans="y"
      endif
   enddo

   open(unit=14,file="freqs.dat")
   write(14,*) 0.0d0
   call freqout(ma,nstar,a,aerr,sns)
   close(14)

endif

xcoo=1.0 !don't try to guess phases in cosinefit.. we know it now!

!make copies of original date, because now we are going to slice and dice
npt2=npt
time2=time
mag2=mag
merr2=merr

allocate(pers2(nstarmax),a2(nfmax),aerr2(nfmax),a3(nfmax),aerr3(nfmax))
pers2=pers
a2=a
aerr2=aerr

sfreq1=1.0/(maxx-minx)
sfreq2=1.0/10.0 !10 days is shortest Period to scan.
slsteps=min(int(ofac*(sfreq2-sfreq1)*npt/nyq),1000000)

open(unit=23,file="frmod.dat")

dnu=(freq2-freq1)/dble(steps)
nu=freq1
ndphase=5
rndphase=dble(ndphase)
allocate(dfreqs(nfmax),dthetas(nfmax),nd(nfmax),dfreqserr(nfmax),       &
 dthetaserr(nfmax))
do i=1,slsteps
   per=1.0d0/nu
   do kk=1,ndphase
      dphase=dble(kk-1)/rndphase
      npt=0 !reset counter
      do j=1,npt2
         phase=time2(j)/per-int(time2(j)/per)+dphase
         if(phase.gt.1.0d0)phase=phase-1.0d0
         if(phase.le.0.5d0)then
!         if((phase.ge.1.0/6.0).and.(phase.lt.2.0/6.0))then
            npt=npt+1
            time(npt)=time2(j)
            mag(npt)=mag2(j)
            merr(npt)=merr(j)
         endif
      enddo
      pers=pers2 !use original solution to seed fit.
      a=a2
      aerr=aerr2
      call pgpanl(1,2)
      call pgeras()
      plot=0 !1 means plot the fit..
!      write(0,*) a(1),pers(1),ma,npt
      call cosinefit(pers,res,kID,xcoo,ycoo,avg,a,btheta,rmavg,plot,    &
        ma,aerr,1,2,w,nstar,fixfreq)
      call fixa(ma,nstar,a)
      a3=a       !save solution
      aerr3=aerr

      npt=0
      do j=1,npt2
         phase=time2(j)/per-int(time2(j)/per)+dphase
         if(phase.gt.1.0d0)phase=phase-1.0d0
         if(phase.gt.0.5d0)then
!         if((phase.ge.4.0/6.0).and.(phase.lt.5.0/6.0))then
            npt=npt+1
            time(npt)=time2(j)
            mag(npt)=mag2(j)
            merr(npt)=merr(j)
         endif
      enddo
      pers=pers2 !use original solution to seed fit.
      a=a2
      aerr=aerr2
      call pgpanl(1,2)
      call pgeras()
      plot=0 !1 means plot the fit..
!      write(0,*) a(1),pers(1),ma,npt
      call cosinefit(pers,res,kID,xcoo,ycoo,avg,a,btheta,rmavg,plot,    &
        ma,aerr,1,2,w,nstar,fixfreq)
      call fixa(ma,nstar,a)

      j=(ma-nstar)/nstar
      do 23 k=1,nstar
         cc=(j+1)*(k-1)+2
         dfreqs(k)=(a3(cc-1)/twopi-a(cc-1)/twopi) !(c/d)
         dfreqserr(k)=sqrt((aerr3(cc-1)/twopi)**2.0+                    &
          (aerr(cc-1)/twopi)**2.0)
!        check both ways of subtracking phase, take min
!         dthetas(k)=min(abs(a3(cc+1)-a(cc+1)),                          &
!          abs(twopi+a3(cc+1)-a(cc+1)),abs(a3(cc+1)-a(cc+1)-twopi))
         dthetas(k)=a3(cc+1)-a(cc+1)
         n=int(dthetas(k)/Pid2)
         dthetas(k)=dthetas(k)-dble(n)*Pid2
         dthetas(k)=dthetas(k)/Pid2/a2(cc-1)*86400.0d0 !(delay is sec)
         dthetaserr(k)=sqrt(aerr3(cc+1)**2.0+aerr(cc+1)**2.0)
         dthetaserr(k)=dthetaserr(k)/Pid2/a2(cc-1)*86400.0d0 !sec
 23   continue

      call rqsort(nstar,dfreqs,nd) !sort data
      jj=nstar/2
      if(jj.le.0) jj=1
      mfreq=dfreqs(nd(jj))  !calculate median
      mstd=mstdev(nstar,dfreqs,mfreq) !calculate standard deviation
      afreq=0.0d0
      abin=0.0d0
      do j=1,nstar
!         write(0,504) dfreqs(j),dfreqserr(j)
         if(abs(dfreqs(j)-mfreq).lt.3.0*mstd)then !cut out outliers
            afreq=afreq+dfreqs(j)/dfreqserr(j)
            abin=abin+1.0/dfreqserr(j)
         endif
      enddo
      afreq=afreq/abin  !weighted mean
!      read(5,*)

      call rqsort(nstar,dthetas,nd)  !sort data
      mtheta=dthetas(nd(jj))  !get median
      mstd=mstdev(nstar,dthetas,mtheta) !calculate std around median
      atheta=0.0d0
      abin=0.0d0
      do j=1,nstar
!         write(0,504) mtheta,mstd,dthetas(j),dthetaserr(j)
         if(abs(dthetas(j)-mtheta).lt.3.0*mstd)then
            atheta=atheta+dthetas(j)/dthetaserr(j)
            abin=abin+1.0/dthetaserr(j)
         endif
      enddo
      atheta=atheta/abin !weighted mean
!      read(5,*)

      write(6,505) "search: ",nu,afreq,atheta
      write(23,504) nu,afreq,atheta  !c/d , c/d, seconds
504   format(17(1PE17.10,1X))
505   format(A8,17(1PE17.10,1X))
   enddo
   nu=nu+dnu
enddo

close(23)

end program freqmod
