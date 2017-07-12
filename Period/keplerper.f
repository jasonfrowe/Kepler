      program keplerper
C     Automatic DFT for Fourier Decomposition
C     Jason Rowe - jasonfrowe@gmail.com
      implicit none
      integer nmax,fixfreq,steps,nfit,iargc,kID,npt,i,nb,nc,bins,nfitl,
     .  nstarmax,niter,nstar,plot,nfmax,ma,itmax,nunit
      parameter(nmax=2000000,nstarmax=300,nfmax=2000)
      integer ndt(nmax)
      double precision Keplertime,tbin,snlimit,freq1,freq2,rmavg,
     .  time(nmax),mag(nmax),merr(nmax),minx,maxx,nyquest,nyq,dt(nmax),
     .  ofac,sig,bper,btheta,sn,pers(nstarmax),w(nstarmax),aerr(nfmax),
     .  a(nfmax),avg,xcoo,ycoo,res(nmax),sns(nstarmax),props(5),boxbin,
     .  std(2),mean,stdev,itime(nmax)
      character*80 cline,ans,title,filename
      logical multi
      common /data/ time,mag,merr,npt,nfit
      
      xcoo=-1.0

      fixfreq=0 !0=fit frequencies, 1=fixed frequencies
      keplertime=2451545.00000000 !HJD of Kepler start time
      tbin=0.0  !initialize binning of data parameters (0==do not bin)
      snlimit=3.6 !default limit for determination of significant amps
      freq1=-1.0 !set to 1/(maxT-minT)
      freq2=-1.0 !set to nyquest
      steps=0 !initialize number of steps (default is ~1200 minimum)
      nfit=1  !for harmonic fits
      
      if(iargc().lt.1) goto 903
      call getarg(1,filename)
      kID=0
c      call getarg(1,cline)
c      read(cline,*) kID  !read in Kepler ID
      if(iargc().ge.2) then
        call getarg(2,cline)
        read(cline,*) freq2 !get high frequency for DFT
      endif
      if(iargc().ge.3) then
         call getarg(3,cline)
         read(cline,*) freq1 !get low frequency for DFT
      endif
      if(iargc().ge.4) then
         call getarg(4,cline)
         read(cline,*) snlimit !get SN-limit for frequency detection
      endif
      if(iargc().ge.5) then
         call getarg(5,cline)
         read(cline,*) tbin !get binning parameter
      endif
      if(iargc().ge.6) then
         call getarg(6,cline)
         read(cline,*) steps !get number of steps of DFT scan
      endif
      if(iargc().ge.7) then
         call getarg(7,cline)
         read(cline,*) nfit !number of harmonics to fit
      endif
      if(iargc().ge.8) then
         call getarg(8,cline)
         read(cline,*) fixfreq !whether or not to fix frequencies
      endif
      if(nfit.lt.1) goto 903
      if((fixfreq.lt.0).or.(fixfreq.gt.2)) goto 904
      multi=.false.
      rmavg=0.0
      
      if(kID.eq.0)then
c        write(6,*) "Enter filename"
c        read(5,*) filename
        nunit=10
        open(unit=nunit,file=filename,status='old',err=901)
        call readkeplc(nunit,nmax,npt,time,mag,merr,itime,
     .  Keplertime)
        close(nunit)
        title=" "
      else
C       Read in Kepler data based on Kepler ID (kID)
        call readdata(kID,nmax,npt,time,mag,merr,Keplertime,props)      
        write(0,*) "kID      Kmag    Teff   logg   rad"
        write(0,500) kID,props(1),props(2),props(3),props(4)
        write(title,500) kID,props(1),props(2),props(3),props(4)
 500    format(I8,1X,F6.3,1X,F7.0,1X,F5.3,1X,F5.2)
      endif
      
c      open(unit=14,file="freqs.dat") !stores frequency fits
      open(unit=22,file="std.dat",POSITION='APPEND')
            
C     These lines convert to parts-per-thousand
      do 11 i=1,npt
        mag(i)=mag(i)*1.0e3
        merr(i)=merr(i)*1.0e3
 11   continue

      write(0,*) "Points Read: ",npt

      minx=time(1)
      maxx=time(1)
      do 5 i=2,npt
         maxx=max(time(i),maxx)
         minx=min(time(i),minx)
 5    continue
 
      call pgopen('?')
c      call pgask(.false.)
      call pgsch(2.9)
      call pgsubp(1,4)
      call pgvport(0.1,0.95,0.2,0.8) !gives enough room for labels

      call pgpage()


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Commands Start Here
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     plotdata
C     pdm
C     plotph
C     cosinefit
C     grelor (not yet functional)
C     windowfn
C     avgrm
C     shaperm

      nyq=nyquest(nmax,npt,time,dt,ndt)

      nb=5
      nc=2
      if(freq1.lt.0.0) freq1=1.0/(maxx-minx)
      if(freq2.lt.0.0) freq2=nyq
c      per1=3.0
c      per2=20.0
c      per1=200.0
c      per2=250.0
      ofac=4.0 !over sampling
c      write(0,*) "nyq",nyq
      if(steps.le.0) steps=int(ofac*(freq2-freq1)*npt/nyq)
c      steps=12000
      if(steps.gt.1e6)then
c         write(0,*) "reducing steps to 1e6"
         steps=1e6
      endif
c      steps=45000
c      steps=4000
c      write(6,*) "freq1,freq2,steps",freq1,freq2,steps
      nfitl=3
      niter=nstarmax-1
      sig=5.0 !set level of sigma clipping

C     sigma clipping
      do 10 i=1,1
        call sigclip(npt,time,mag,merr,sig,rmavg)
 10   continue
 
      bins=150 !number of bins for primative binning
c      call bind(npt,time,mag,merr,bins) !primative binning (ugh!)
      call avgrm(npt,mag,rmavg) !remove zero-point from data.
      if(tbin.gt.0.0) call bindt(npt,time,mag,merr,tbin) !not tested!!
c      call detrend(npt,time,mag,merr,nfitl) !polynomial detrending
      boxbin=1.0 !boxcar filter in days
c      call boxcar(npt,time,mag,merr,boxbin)
c      tbin=1.0*500.0
      tbin=0.25*24.0*60.0
c      call spdetrend(npt,time,mag,merr,tbin,sky) !spline detrending
c      call sigclip(npt,time,mag,merr,sig,rmavg)
c      goto 21
c      call plotdata(npt,time,mag,merr,-1.0,-1.0,1,1,MOSTtime)
c      write(6,*) id,xcoo,ycoo,npt
      call avgrm(npt,mag,rmavg) !make sure zero point is up-to-date
      
      std(1)=1.0d3*stdev(npt,mag,mean) !standard deviation in ppm

c      write(0,*) "rmavg: ",rmavg

      nstar=0
      
      write(0,*) freq1,freq2,steps
      
c      call plotph(npt,time,mag,3.52474859,1,2)
      call plotdata(npt,time,mag,merr,-1.0d0,-1.0d0,1,1,Keplertime,
     .  title)
      call jmfourw(npt,time,mag,merr,freq1,freq2,steps,bper,btheta,sn,
     .  1,2,snlimit,1)
c      goto 31

      write(0,*) "bper:",bper
      plot=0
      bper=1.0/bper !cosinefit works with Periods not frequencies 
      nstar=nstar+1 !number of individual frequencies for fitting.
      pers(nstar)=bper !store best *period*
      sns(nstar)=sn !store SN of best period
      call cosinefit(pers,res,kid,xcoo,ycoo,avg,a,btheta,rmavg,plot,
     .     ma,aerr,1,3,w,nstar,fixfreq)
      call windowfn(npt,nfit,time,mag,merr,bper,freq1,freq2,avg,a,nb,nc,
     .   steps,1,3,nstar,pers)
c      goto 31
      plot=0
      call shaperm(npt,nfit,time,mag,merr,res,pers,avg,a,1,3,plot,w,
     .     nstar,Keplertime)
      call jmfourw(npt,time,res,merr,freq1,freq2,steps,bper,btheta,sn,
     .  1,4,snlimit,1)
     
      plot=1
      write(0,*) "Next peak?          "
      read(5,501) ans
 501  format(A1)
      if(ans.eq."n")goto 24
c      ans="y"

      itmax=50
      
      i=0

  20  i=i+1
      if((sn.lt.snlimit).or.(i.gt.itmax)) then
         ans="n"
      else
         ans="y"
      endif

      if(ans.eq."y")then
         write(0,*) "Iteration #: ",i
         write(0,*) "bfreq:",bper,sn
cp         call pgpanl(1,4)
cp         call pgpage()
         bper=1.0/bper
         nstar=nstar+1
         pers(nstar)=bper
         sns(nstar)=sn
cp         call jmfour(npt,time,mag,per1,per2,steps,bper,btheta,1,1)
         call pgpanl(1,2)
         call pgeras()
c         call plotph(npt,time,res,bper,1,2)
         plot=1
         call cosinefit(pers,res,kID,xcoo,ycoo,avg,a,btheta,rmavg,plot,
     .        ma,aerr,1,2,w,nstar,fixfreq)
         write(0,503) i+1,"bfreq:",1.0/bper,sn
 503     format(i3,1x,A6,1X,F8.5,1X,F6.3)
c         call plotdata(npt,time,mag,merr,-1.0,-1.0,1,2,MOSTtime)
c         call plotfit(npt,nfit,time,mag,pers,avg,a,w,nstar)
         call pgpanl(1,3)
         call pgeras()
         call shaperm(npt,nfit,time,mag,merr,res,pers,avg,a,1,3,plot,w,
     .        nstar,Keplertime)
         call pgpanl(1,4)
c         call pgeras()
         call jmfourw(npt,time,res,merr,freq1,freq2,steps,bper,btheta,
     .      sn,1,4,snlimit,1)
c         call pdm(npt,time,res,nb,nc,per1,per2,steps,bper,btheta,1,4)
         goto 20
      endif

 
  24  continue

c      call pgpanl(1,1)
c      call pgeras()
c      call plotdata(npt,time,mag,merr,-1.0d0,-1.0d0,1,1,Keplertime,
c     .  title)
c      call plotfit(npt,nfit,time,mag,pers,avg,a,w,nstar)
      
      open(unit=14,file="freqs.dat")
      write(14,*) 0.0d0
      call freqout(ma,nstar,a,aerr,sns)
      close(14)

C     Standard deviation of pre-whitened data in ppm
      std(2)=1.0d3*stdev(npt,res,mean)
      
      do 30 i=1,2
        if(std(i).ge.100000.0) std(i)=99999.0
 30   continue

      
      write(22,502) kID,props(1),props(2),props(3),props(4),props(5),
     .  std(1),std(2),nstar
 502  format(I8,1X,F6.3,1X,F7.0,1X,F5.3,1X,F5.2,1X,F3.0,2(1X,F8.2),
     .  1X,I3)
      close(22)
      
      open(unit=23,file='test.dat')
      do 32 i=1,npt
        write(23,*) time(i)+Keplertime,mag(i),merr(i)
 32   continue
      close(23)
     
C     close(14)

 31   call pgclos()
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 903  write(0,*) "Usage: keplerper <KeplerID> <freq2> [freq1] [sn] [tbin
     .] [steps] [nfit] [fixfreq]"
      write(0,*) "<KeplerID>: KeplerID of photometry to extract"
      write(0,*) "<freq2>: largest frequency to scan"
      write(0,*) "[freq1]: lowest frequency to scan (optional)"
      write(0,*) "         set freq1 negative for default"
      write(0,*) "[sn]: S/N for statistically significant frequencies (o
     .ptional)"
      write(0,*) "        default = 3.6"
      write(0,*) "[tbin]: Binning time scale (minutes) (optional)"
      write(0,*) "        set to < 0 to disable binning"
      write(0,*) "[steps]: Steps for DFT scan"
      write(0,*) "         default is 4 times oversampling"
      write(0,*) "[nfit]: number of harmonics to fit (optional)"
      write(0,*) "        default is 1"
      write(0,*) "[fixfreq]: 0 (default): fit frequencies"
      write(0,*) "           1 : do not fit frequencies"
      write(0,*) "           2 : do not fix frequencies, except for fina
     .l iteration"
      goto 999
 904  write(0,*) "Bad Usage:"
      write(0,*) "[fixfreq] must be 0, 1 or 2."
      goto 999  
 999  end
