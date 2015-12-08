      program transitfinder
      implicit none
      integer iargc,nunit,nmax,npt,i,steps,nb,in1,in2,nbf,i1,i2,nplot,j,
     .  k,nsamp
      parameter(nmax=600000)
      integer ndt(nmax),nd(nmax)
      double precision time(nmax),flux(nmax),ferr(nmax),itime(nmax),
     .  MOSTtime,u(nmax),v(nmax),nyq,nyquest,dt(nmax),maxx,minx,ofac,
     .  freq1,freq2,df,qmi,qma,p(nmax),bper,bpow,depth,qtran,epo,pmean,
     .  bfreq,freqs(nmax),median,std,stdev,pcut(nmax),perr(nmax),
     .  psmooth(nmax),phase,sn(nmax),fread,phase_off,phase_start,
     .  phase_end,tphase
      character*80 obsfile,cline,txtout
      
      
      if(iargc().lt.1) goto 901 !check number of command line arguements
      call getarg(1,obsfile) !get filename for input data
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=901)
c      call readdata(nunit,nmax,npt,time,flux,ferr,itime,MOSTtime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,MOSTtime)
      close(nunit)!release unit number as we are done with the file.
      
      nplot=0 !disable plotting by default
      if(iargc().ge.2)then
        call getarg(2,cline)
        read(cline,*) nplot
      endif
      
      minx=time(1)
      maxx=time(1)
      do 10 i=2,npt
        minx=min(minx,time(i))
        maxx=max(maxx,time(i))
 10   continue


      nyq=nyquest(nmax,npt,time,dt,ndt)
c      write(0,*) "Nyquest: ",nyq
      
      freq1=2.0/(maxx-minx)
      if(iargc().ge.3)then
        call getarg(3,cline)
        read(cline,*) fread
        if(fread.gt.0.0d0) freq1=fread
      endif

c      freq2=2.0 !0.5 day min period
      freq2=(1.0/2.0)+0.01 !3.0 day min period
      freq2=0.45
c      write(0,*) "F:",freq1,freq2

      if(iargc().ge.4)then
        call getarg(4,cline)
        read(cline,*) freq2
      endif


      ofac=1024.0 !over sampling
      
      steps=int(ofac*(freq2-freq1)*npt/nyq)
c      write(0,*) "steps:",steps
      df=(freq2-freq1)/dble(steps)

      write(0,*) "Steps: ",steps
      nsamp=max(int(ofac),steps/1000) !sample size of std,median.. etc

      nb=1900
      qmi=0.00030
c      qma=0.05
      qma=0.005
      
      do 19 i=1,steps
        sn(i)=1.0d0
        psmooth(i)=0.0d0
 19   continue

      call eebls(npt,time,flux,u,v,steps,freq1,df,nb,qmi,qma,p,bper,
     .  bpow,depth,qtran,in1,in2,sn,psmooth)
     
      do 12 i=1,steps
        freqs(i)=freq1+(i-1)*df
        perr(i)=0.001d0
 12   continue
      
      do 13 i=1,steps
        i1=max(1,i-nsamp)
        i2=min(steps,i+nsamp)
        k=0
        do 14 j=i1,i2
            k=k+1
            pcut(k)=p(j)
 14     continue
        psmooth(i)=pcut(1)
        do 16 j=2,k
            psmooth(i)=min(psmooth(i),pcut(j))
 16     continue
c        call rqsort(k,pcut,nd)
c        psmooth(i)=pcut(nd(k/2))
 13   continue
      
      do 15 i=1,steps
        p(i)=p(i)-psmooth(i)
!c------------------
!csn        i1=max(1,i-nsamp)
!csn        i2=min(steps,i+nsamp)
!c        sn(i)=0.0
!c        do 18 j=i1,i2
!c            sn(i)=sn(i)+p(j)
!c 18     continue
!c        sn(i)=sn(i)/dble(i2-i1+1)
!csn        k=0
!csn        do 21 j=i1,i2
!csn            k=k+1
!csn            pcut(k)=p(j)
!csn 21     continue
!csn        call rqsort(k,pcut,nd)
        sn(i)=1.0!stdev(j,pcut,pmean)!1.0d0!pcut(nd(k/2))
!c------------------
 15   continue
 
C     ----ADD THIS FOR FREQUENCY SMOOTHING
      call eebls(npt,time,flux,u,v,steps,freq1,df,nb,qmi,qma,p,bper,
     .  bpow,depth,qtran,in1,in2,sn,psmooth)
!      do 20 i=1,steps
!        p(i)=p(i)-psmooth(i)
! 20   continue
 
!      p(1)=0.0d0
     
!      bpow=p(1)
!      bper=1.0d0/freq1
!      do 17 i=2,steps-1
!        p(i)=p(i)/sn(i)
!        if(p(i).gt.bpow)then
!            bpow=p(i)
!            bper=1.0d0/(freq1+(i-1)*df)
!        endif
! 17   continue
     
!      write(0,*) "in#",in1,in2,minx
      bfreq=1.0/bper
      nbf=(bfreq-freq1)/df +1
C     Time of Transit
      if(in1.lt.in2)then
        epo=minx+bper*dble((in1+in2)/2.0)/dble(nb)
      else
        epo=minx+bper*dble((in1+in2-nb)/2.0)/dble(nb)
      endif
      
c      i1=max(1,in2+1)
c      i2=min(steps,in2+int((in2-in1+1)*(qtran+1)/qtran))
c      write(0,*) i1,i2

c Calculating the Standard Deviation
c PHASE OFFSET
      phase_off = (epo/bper) - int(epo/bper)
c      write(0,*) phase_off
c START AND END PHASE
      phase_start = 1.0 - (2.0*qtran) - int(2.0*qtran)
      phase_end = (2.0*qtran) - int(2.0*qtran)
      if(bper*qtran .gt. 0.5*bper)then
         phase_start = 0.75
         phase_end = 0.25
      endif
c only want data outside the transit
      j = 0 !initialize
      do 11 i = 1, npt
      tphase = (time(i)/bper) - int(time(i)/bper) - phase_off
      if(tphase .lt. 0)then
         tphase = tphase + 1.0
      endif
      if((tphase.gt.phase_end).and.(tphase.lt.phase_start))then
         j = j + 1
         pcut(j)=flux(i)
      endif
11    continue
c      write(0,*) j
      std=stdev(j, pcut, pmean)
      call rqsort(j,pcut,nd)
      median=pcut(nd(j/2))
c call standard deviation routine

c      std=stdev(npt,flux,pmean)
c      call rqsort(npt,flux,nd)
c      median=flux(nd(npt/2))

c      write(0,*) i1,i2,steps,std,median
      
c      write(6,*) minx,maxx
c      write(0,*) qtran,std,depth
      write(6,500) bper,epo,bpow,
     .   (depth-pmean)/std*sqrt(qtran*dble(npt)),qtran*bper
 500  format(6(1PE17.10,1X))
c      write(6,*) in1,in2

      if(nplot.ge.1)then
        if(nplot.eq.1) call pgopen('/xserve')
        if(nplot.eq.2) call pgopen('find.ps/vcps')
        if(nplot.eq.3) call pgopen('find.png/png')
c       call pgask(.false.)
        call pgsch(2.9)
        call pgsubp(1,4)
        call pgvport(0.1,0.95,0.2,0.8) !gives enough room for labels

        call pgpage()
        
        call plot(steps,freqs,p,obsfile,"frequency c/d","P",0,0.0,0.0)
        call pgpage()
        call plot(npt,time,flux,obsfile,"BJD-2454900","Flux",1,epo,bper)
        call pgpage()
        phase=epo/bper-int(epo/bper)
        if(phase.lt.0.0d0) phase=phase+1.0d0
        call plotph(npt,time,flux,bper,phase)
        call pgpage()
        call pgwindow(0.0,1.0,0.0,1.0)

        write(txtout,501) "Per: ",bper
 501    format(A5,1X,1PE17.10)
        call pgptxt(0.1,0.9,0.0,0.5,txtout)
        write(txtout,501) "Epo: ",epo
        call pgptxt(0.1,0.77,0.0,0.5,txtout)
        write(txtout,501) "Pow: ",bpow
        call pgptxt(0.1,0.64,0.0,0.5,txtout)
       write(txtout,501) "S/N: ",(depth-pmean)/std*sqrt(qtran*dble(npt))
        call pgptxt(0.1,0.51,0.0,0.5,txtout)
        write(txtout,501) "Dur: ",qtran*bper
        call pgptxt(0.1,0.38,0.0,0.5,txtout)
        
        call plottrans(npt,time,flux,bper,phase,qtran)

        call pgclos()
      endif
      
      goto 999
 901  write(0,*) "Usage: transitfind <photfile> [np] [freqlow] [freqhi]"
      write(0,*) " photfile : photometry file"
      write(0,*) " np : 0-no plot, 1-xwindow, 2-ps file, 3-png file"
      write(0,*) " freqlow,freqhi : frequency range to scan"
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plottrans(n,x,y,per,phase,qtran)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,n
      real bb(4),px,py
      double precision x(n),y(n),per,phase,qtran

      call pgvport(0.5,0.95,0.2,0.8) !gives enough room for labels

      bb(1)=real(max(0.5-2.0*qtran,0.0))
      bb(2)=real(min(0.5+2.0*qtran,1.0))
      bb(3)=y(1)
      bb(4)=y(1)
      do 10 i=1,n
        bb(3)=min(bb(3),y(i))
        bb(4)=max(bb(4),y(i))
 10   continue

      call pgwindow(bb(1),bb(2),bb(3),bb(4))
c      call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)
      call pgptxt(bb(1)-0.15*(bb(2)-bb(1)),(bb(3)+bb(4))/2.0,90.0,0.5,
     .   "flux")
      call pgptxt((bb(2)+bb(1))/2.0,bb(3)-0.25*(bb(4)-bb(3)), 0.0,0.5,
     .   "phase (hours)")

      call pgbbuf()
      do 30 i=1,n
        px=real(x(i)/per-int(x(i)/per)-phase+0.5)
        if(px.lt.0.0)px=px+1
        py=real(y(i))
        call pgpt1(px,py,-1)
 30   continue
      call pgebuf()

      call pgwindow(real(-48.0*qtran*per),real(48.0*qtran*per),
     .   bb(3),bb(4))
      call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotph(n,x,y,per,phase)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n,i,stepmax
      parameter(stepmax=1200)
      real px,py,bb(4),lx(2),ly(2)
      double precision x(n),y(n),per,phase
      
      bb(1)=0.0
      bb(2)=1.0
      bb(3)=y(1)
      bb(4)=y(1)
      do 10 i=1,n
        bb(3)=min(bb(3),y(i))
        bb(4)=max(bb(4),y(i))
 10   continue
   
      call pgwindow(bb(1),bb(2),bb(3),bb(4))
      call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)
      call pgptxt(bb(1)-0.08*(bb(2)-bb(1)),(bb(3)+bb(4))/2.0,90.0,0.5,
     .   "flux")
      call pgptxt((bb(2)+bb(1))/2.0,bb(3)-0.25*(bb(4)-bb(3)), 0.0,0.5,
     .   "phase")
      
      lx(1)=real(phase)
      lx(2)=real(phase)
      ly(1)=bb(3)
      ly(2)=bb(4)
      call pgsci(2)
      call pgline(2,lx,ly)
      call pgsci(1)
      
      call pgbbuf()
      do 30 i=1,n
        px=real(x(i)/per-int(x(i)/per))
        if(px.lt.0.0)px=px+1
        py=real(y(i))
        call pgpt1(px,py,-1)
 30   continue
      call pgebuf()
 
      return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plot(n,x,y,title,xlabel,ylabel,nmark,epo,per)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n,i,stepmax,nmark
      parameter(stepmax=1200)
      real px(3),py(3),bb(4),xp,yp
      double precision x(n),y(n),epo,per
      character*(*) title,xlabel,ylabel
      
      bb(1)=x(1)
      bb(2)=x(1)
      bb(3)=y(1)
      bb(4)=y(1)
      do 10 i=1,n
        bb(1)=min(bb(1),x(i))
        bb(2)=max(bb(2),x(i))
        bb(3)=min(bb(3),y(i))
        bb(4)=max(bb(4),y(i))
 10   continue
   
      call pgwindow(bb(1),bb(2),bb(3),bb(4))
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel(xlabel,ylabel,title)
      call pgptxt(bb(1)-0.08*(bb(2)-bb(1)),(bb(3)+bb(4))/2.0,90.0,0.5,
     .   ylabel)
      call pgptxt((bb(2)+bb(1))/2.0,bb(3)-0.23*(bb(4)-bb(3)), 0.0,0.5,
     .   xlabel)
      call pgptxt((bb(2)+bb(1))/2.0,bb(4)+0.10*(bb(4)-bb(3)), 0.0,0.5,
     .   title)
      
      call pgbbuf()
      do 30 i=2,n
        px(1)=real(x(i-1))
        py(1)=real(y(i-1))
        px(2)=real(x(i))
        px(3)=real(x(i))
        py(2)=real(y(i))
        py(3)=real(y(i))
        call pgline(3,px,py)
 30   continue
      call pgebuf()
 
      if(nmark.eq.1)then
         xp=bb(1)
         i=0
         yp=bb(3)+0.05*(bb(4)-bb(3))
         do while((xp.lt.bb(2)).or.(i.lt.10000))
            xp=real(epo+per*dble(i))
            call pgsci(2)
            call pgpt1(xp,yp,13)
            call pgsci(1)
            i=i+1
         enddo
      endif

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function stdev(npt,pts,mean)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculates standard deviation of data set given the mean.
      implicit none

      integer npt,i
      double precision pts(npt),mean,s,ep,adev,p,var,sdev,mean2

      s=0.
      do 11 i=1,npt
         s=s+pts(i)
 11   continue
      mean=s/npt

      ep=0.
      var=0.
      do 10 i=1,npt
         s=pts(i)-mean
         ep=ep+s
         p=s*s
         var=var+p
 10   continue
      var=(var-ep**2/npt)/(npt-1)
      sdev=sqrt(var)

      stdev=sdev

      return
      end
      

        
