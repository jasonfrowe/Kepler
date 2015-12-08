CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine jmfourw(npt,time,mag,merr,freq1,freq2,steps,bper,
     .  btheta,signoise,panx,pany,snlimit,plot)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,steps,panx,pany,nmax,stepmax,k,j,nf,snw,n1,n2,plot,
     .   ne6
      parameter(stepmax=2000,snw=2000,ne6=1000000)
      parameter(nmax=2000000)
      real percold,perc,bb(4),px(3),py(3),orbper,px2(ne6),py2(ne6),dumr
      double precision time(npt),mag(npt),freq1,freq2,bper,btheta,by(4),
     .  nx(stepmax),ny1(stepmax),ny2(stepmax),bx(4),tbin,dum,weight,
     .  merr(npt),amps(ne6),cf(2*snw+2),stdev,std,cf2(2*snw+2),mean,
     .  signoise,snplot(ne6),freqs(ne6),snlimit
      character*80 tmpc

      INTEGER I,NUMPTS,NMDENU
      DOUBLE PRECISION JD (nmax),B(nmax),LONU,DELTNU,TWOPI
      DOUBLE PRECISION SIDEL(nmax),CODEL(nmax),SI(nmax),CO(nmax)
      DOUBLE PRECISION SISUM,COSUM,SN,NU,AMP,PHI
c      CHARACTER *20 INFILE,SPFILE
      DOUBLE PRECISION JDTMP,LOJD,HIJD,BTMP,JDZRO

c      snlimit=3.6

      if(steps.lt.stepmax) then
c         write(6,*) "Increasing number of steps to ",stepmax
         steps=stepmax
      endif


      tbin=0.0

      btheta=-99.9d10
      do 5 i=1,stepmax
         ny1(i)=-99.9d10
         ny2(i)= 99.9d10
 5    continue
C     Uncomment the c21 lines to output the FT to a text file.
c      open(unit=21,file="ft.dat")

      LOJD=0.
      HIJD=100.
      JDZRO=0.

      LONU=freq1
      DELTNU=(freq2-freq1)/dble(steps)
      NMDENU=steps

      do 3 i=1,npt
         jd(i)=time(i)
         b(i)=mag(i)
 3    continue
      NUMPTS = npt

      TWOPI = 6.2831853
      NU = LONU
      DO 61 I=1,NUMPTS
        SI(I) = SIN(TWOPI*LONU*JD(I))
        CO(I) = COS(TWOPI*LONU*JD(I))
        SIDEL(I) = SIN(TWOPI*DELTNU*JD(I))
        CODEL(I) = COS(TWOPI*DELTNU*JD(I))
   61 CONTINUE

      percold=0.0      
      DO 51 J=1,NMDENU
c         perc=100.0*real(j)/real(steps)
c         if(perc-percold.gt.5.0) then
c            write(tmpc,500) "Percent done ", perc
c            percold=perc
c            call ovrwrt(tmpc,2)
c         endif
 500     format(A13,F6.1)
         SISUM = 0.
         COSUM = 0.
         weight=0.0
         DO 41 I=1,NUMPTS
            SISUM = SISUM + B(I)*SI(I)/merr(i)
            COSUM = COSUM + B(I)*CO(I)/merr(i)
            weight=weight+1.0/merr(i)
 41      CONTINUE
C     WRITE (50,*) NU,SISUM,COSUM
         AMP = (2./weight)*SQRT(SISUM**2. + COSUM**2.)
         PHI = ATAN(-SISUM/COSUM)
C     Doh.. mistake
c         PHI= ATAN(COSUM/SISUM)
         IF (COSUM.LT.0) PHI = PHI + TWOPI/2.
         IF ((SISUM.GT.0).AND.(COSUM.GE.0)) PHI = PHI + TWOPI
         amps(j)=amp
         freqs(j)=nu
         k=1+dble(stepmax)/dble(steps)*(j-1)
         nx(k)=nu
         if(amp.gt.ny1(k)) ny1(k)=amp
         if(amp.lt.ny2(k)) ny2(k)=amp
c         if(tbin.gt.0.) then
c         if((amp.gt.btheta).and.(nu.gt.1440.0/tbin)) then
c            btheta=amp
c            bper=nu
c            nf=j
c         endif
c         else
         if(amp.gt.btheta) then
            btheta=amp
            bper=nu
            nf=j
         endif
c         endif 
c     call pgpt1(NU,AMP,1)
c21      WRITE (21,*) NU,AMP,PHI
c     
         DO 31 I=1,NUMPTS
            SN = SI(I)
            SI(I) = SN*CODEL(I) + CO(I)*SIDEL(I)
            CO(I) = CO(I)*CODEL(I) - SN*SIDEL(I)
 31      CONTINUE
         NU = NU + DELTNU
 51   CONTINUE
 
cC     do s/n calc.
c      n1=nf-snw-1
c      n2=nf+snw
c      if(n1.le.0)n1=1
c      if(n2.gt.NMDENU)n2=NMDENU
c      j=0
c      do 200 i=n1,n2
c         j=j+1
c         cf(j)=amps(i)
c 200  continue
c      std=stdev(j,cf,mean)
c      k=0
c      do 201 i=1,j
c         if(abs(cf(i)-mean).lt.3.0*std) then
c            k=k+1
c            cf2(k)=cf(i)
c         endif
c 201  continue
c      std=stdev(k,cf2,mean)
c      signoise=btheta/mean
      
      bb(1)=real(freq1)
      bb(2)=real(freq2)
      bb(3)=0.0
      bb(4)=real(btheta)
      
      if(plot.eq.1)then
c      call pgpage()
        call pgpanl(panx,pany)
        call pgeras()
        call pgwindow(bb(1),bb(2),bb(3),bb(4))
cc      call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pgbox('BTS',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel("Frequency","AMP"," ")
        call pgptxt((bb(1)+bb(2))/2.0,-0.25*(bb(4)),0.0,0.5,
     .      "Frequency (c/d)")
        call pgptxt((bb(1)+bb(2))/2.0,bb(4)+0.23*(bb(4)),0.0,0.5,
     .      "Frequency (\(0638)Hz)")
        call pgptxt(bb(1)-0.05*(bb(2)-bb(1)),(bb(4))/2,90.0,0.5,
     .      "Amplitude (ppm)")
        call pgbbuf()
      
c21      close(21)
        do 30 i=2,stepmax
            px(1)=real(nx(i-1))
            py(1)=real(ny2(i-1))
            px(2)=real(nx(i))
            px(3)=real(nx(i))
            py(2)=real(ny1(i))
            py(3)=real(ny2(i))
            call pgline(3,px,py)
 30     continue
        call pgebuf()
      endif


C     ***Plotting orbital period and harmonics
c      call pgsls(2)
c      call pgsci(2)
c      orbper=14.199363
c      i=bb(1)/orbper+1
c      j=bb(2)/orbper
c      
c      do 100 k=i,j
c         px(1)=k*orbper
c         px(2)=px(1)
c         py(1)=0.0
c         py(2)=real(btheta)
c         call pgline(2,px,py)
c 100  continue
c
c      call pgsls(1)
c      call pgsci(1)
      
      signoise=0.0
      btheta=0.0
C     do s/n plot
      do 302 nf=1,steps
         n1=nf-snw-1
         n2=nf+snw
         if(n1.le.0)n1=1
         if(n2.gt.steps)n2=steps
         j=0
         do 300 i=n1,n2
            j=j+1
            cf(j)=amps(i)
 300     continue
         std=stdev(j,cf,mean)
         k=0
         do 301 i=1,j
            if(abs(cf(i)-mean).lt.3.0*std) then
               k=k+1
               cf2(k)=cf(i)
            endif
 301     continue
         std=stdev(k,cf2,mean)
         snplot(nf)=snlimit*mean
         if(amps(nf)/mean.gt.snlimit)then  !is S/N good?
c            if((amps(nf).gt.btheta).and.(freqs(nf).le.1.0d0))then  !pick largest amplitude
            if(amps(nf).gt.btheta)then
               signoise=amps(nf)/mean
               bper=freqs(nf)
               btheta=amps(nf)
            endif
         endif
c         write(6,*) n1,n2,nf,freqs(nf),mean
c         read(5,*)
 302  continue
 
      if(plot.eq.1)then
        call pgsci(3)
c      call pgsls(2)
C     draw S/N line
        do 303 i=1,steps
            px2(i)=real(freqs(i))
            py2(i)=real(snplot(i))
 303    continue
        call pgline(steps,px2,py2)
        call pgsci(1)
c      call pgsls(1)

        call pgwindow(bb(1),bb(2),bb(3),1000.0*bb(4))
c        call pgwindow(bb(1),bb(2),bb(3),bb(4))
        call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)

        dumr=1.0e6/86400.0
        call pgwindow(dumr*bb(1),dumr*bb(2),bb(3),1000.0*bb(4))
c        call pgwindow(dumr*bb(1),dumr*bb(2),bb(3),bb(4))
        call pgbox('CMTS1',0.0,0,'',0.0,0)
      endif
c      call pgbox('CTS',0.0,0,'',0.0,0)

c      px(1)=360.0/0.997268
c      px(2)=px(1)
c      py(1)=0.0
c      py(2)=btheta
c      call pgline(2,px,py)
c      call pgsls(1)
c      if(tbin.gt.0.0) then
c         bx(1)=per1
c         by(1)=0.0
c         bx(2)=(60.0*24.0)/tbin
c         by(2)=0.0
c         bx(3)=bx(2)
c         by(3)=btheta
c         bx(4)=per1
c         by(4)=btheta
c         call pgsfs(3)
c         call pgpoly(4,bx,by)
c         call pgsfs(1)
c      endif
c      call pgsci(1)
      
      return
      end
