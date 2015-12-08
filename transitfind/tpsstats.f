      program tpsclean
      implicit none
      integer nmax,nunit,i,j,flag,nunit2,npt,nunit3,nkoi,nEB,nVar,
     .  nunit4,nKIC,nunit5,nbin,nbinmax,iargc
      parameter(nmax=3000000,nbinmax=500)
      integer np(nmax)
      integer kid(nmax),koi(nmax),eb(nmax),var(nmax),kic(nmax)
      real rp(nmax),bdatax(nbinmax),bdatay(nbinmax)
      double precision mes(nmax),ses(nmax),per(nmax),epo(nmax),tce(5),
     .  dumr,kmag(nmax),mag(nmax),x(nmax),y(nmax),dp2(nmax)
      character*80 filename,filenameEB,filenameKID,filenameVAR,
     .  filenameKIC,title
      
      nunit=10
c      filename="tps_q1q3.dat"
      if(iargc().lt.1) goto 907
      call getarg(1,filename) 

      nunit2=11
      filenameEB="eb.list"
      nunit3=12
      filenameKID="KOI.list"
      nunit4=13
      filenameVAR="var.list"
      nunit5=14
      filenameKIC="kic_colour_16.5.txt"
      
      
      open(unit=nunit3,file=filenameKID,status='old',err=904)      
      i=1
 19   read(nunit3,*,end=20) dumr,koi(i)
        i=i+1
      goto 19
 20   continue
      nkoi=i-1
      close(nunit2)

      open(unit=nunit2,file=filenameEB,status='old',err=902)
      i=1
 24   read(nunit2,*,end=23) eb(i)
        i=i+1
      goto 24
 23   continue
      nEB=i-1
      close(nunit2)
      
      open(unit=nunit4,file=filenameVAR,status='old',err=905)
      i=1
 29   read(nunit4,*,end=28) var(i)
        i=i+1
      goto 29
 28   continue
      nVAR=i-1
      close(nunit4)

      open(unit=nunit5,file=filenameKIC,status='old',err=906)
      i=1
 30   read(nunit5,*,end=31) kic(i),kmag(i)
        i=i+1
      goto 30
 31   continue
      nKIC=i-1
      close(nunit5)
      
      open(unit=nunit,file=filename,status='old',err=901)
      
      i=1      
      read(nunit,*) (tce(j),j=1,5)
      kid(i)=int(tce(1))
      per(i)=tce(2)
      ses(i)=tce(3)
      mes(i)=tce(4)
      epo(i)=tce(5)
      
 11   read(nunit,*,end=12) (tce(j),j=1,5)
        flag=0
        if(tce(4).lt.7.0) flag=1
        if(tce(3).lt.sqrt(2.0)) flag=1
        if(flag.eq.0)then
            do 10 j=1,i  !scan through existing KOIs
                if(kid(j).eq.int(tce(1)))then !does KOI already exist?
                    flag=1
                    if(tce(4).gt.mes(j))then !if better, take it
                        kid(j)=int(tce(1))
                        per(j)=tce(2)
                        ses(j)=tce(3)
                        mes(j)=tce(4)
                        epo(j)=tce(5)
                    endif
                endif
 10         continue
        endif
        if(flag.eq.0)then !if it's a new KOI...
            i=i+1
            kid(i)=int(tce(1))
            per(i)=tce(2)
            ses(i)=tce(3)
            mes(i)=tce(4)
            epo(i)=tce(5)
c            write(6,500) kid(i),per(i),ses(i),mes(i),epo(i)
        endif
      goto 11
 12   continue
      npt=i
      
      j=0
      do 25 i=1,npt
        flag=0
c        do 26 j=1,nEB
c            if(eb(j).eq.kid(i)) flag=1
c 26     continue
c        do 27 j=1,nKOI
c            if(koi(j).eq.kid(i)) flag=1
c 27     continue
c        do 32 j=1,nVAR
c            if(var(j).eq.kid(i)) flag=1
c 32     continue
c        do 33 j=1,nKIC
c            if(KIC(j).eq.kid(i)) mag(i)=kmag(j)
c 33     continue
        if(flag.eq.0) then
c            write(6,500) kid(i),per(i),ses(i),mes(i),epo(i),mag(i)
            j=j+1
            x(j)=per(i)
            y(j)=mes(i)
        endif
 25   continue
 500  format(I8,6(1X,F12.6))
 
      call pgopen('?') !open PGPlot device
      call pgask(.true.) !don't ask for new page.. just do it.
      call PGPAP ( 6.0 ,1.0) !paper size
c      call pgsubp(5,5)  !break up plot into grid
      
      nbin=30
      title='Period (days)'
      call histogram(j,rp,per,dp2,np,nbin,nbinmax,bdatax,bdatay,title)
      
      call pgclos()
      
      close(nunit)
      
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 902  write(0,*) "Cannot open ",filenameEB
      goto 999     
 903  write(0,*) "error on line",j
      goto 999
 904  write(0,*) "Cannot open ",filenameKID
      goto 999
 905  write(0,*) "Cannot open ",filenameVAR
      goto 999
 906  write(0,*) "Cannot open ",filenameKIC
      goto 999   
 907  write(0,*) "Usage: tpsstats tpsoutput"
      goto 999 
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine histogram(npt,rp,dp,dp2,np,nbin,nbinmax,bdatax,bdatay,
     .  title)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,j,npt2,nbin,nbinmax
      integer np(npt)
      real rp(npt),datamin,datamax,bdatax(nbinmax),bdatay(nbinmax),bmax,
     .  rave,errs(6),std,rmed
      double precision dp(npt),dp2(npt),ave,var,med
      character*80 title
      
C     First we remove the average (helps a lot with double -> real)
      call avevar(dp,npt,ave,var)
      rave=real(ave) !convert to real*4
      std=real(sqrt(var))
      
C     Now we convert dp to rp
      j=0 !because we have sigma-clipping, we need a counter.
      do 10 i=1,npt
c        if(abs(dp(i)-ave).lt.4.0*sqrt(var))then
            j=j+1
            dp2(j)=dp(i)-ave
            rp(j)=real(dp2(j))
c        endif
 10   continue
      npt2=j

cC     Find median          
c      call rqsort(npt2,dp2,np) !changed npt to k
c      i=npt2/2
c      if(i.le.0) i=1
c      med=dp2(np(i))
c      rmed=real(med)

C     Find datarange
      datamin=rp(1)
      datamax=rp(1)
      do 12 i=2,npt2
        datamin=min(rp(i),datamin)
        datamax=max(rp(i),datamax)
 12   continue
 
      call bindata(nbin,npt2,rp,bdatax,bdatay,datamin,datamax,bmax)
      
      call pgpage() !fresh plotting surface
      call pgslw(1)
      call pgsch(2.0)
c         call windowsetup(xb1,xb2,yb1,yb2) !make a square plotting surface
c         call pgvport(xb1,xb2,yb1,yb2)
      call pgvport(0.2,1.0,0.2,0.9)
c      fc=2.*(datamax-datamin)/real(nbin) !center histograms
      call pgwindow(datamin,datamax,0.,bmax+0.1*bmax) !set size
         
C        Add axis labels
      call pglabel(title,"Relative Probability","")
      call pgbin(nbin,bdatax,bdatay,.false.) !plot the histogram
      
C     Shift axis scale to account for average removal
      call pgwindow(datamin+rave,datamax+rave,0.,1.0+0.1*1.0)
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0) !add boarders
      
c      call errorest(npt2,rp,nbin,bdatax,bdatay,bmax,rmed,errs)
cc      call errorest2(nbin,bdatax,bdatay,bmax,rmed,errs)
c      rmed=rmed+rave !correct for average removal
cC     Need to recalulate standard deviation at this point!      

      return
      end      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindata(nbin,npt,pdata,bdatax,bdatay,datamin,datamax,
     .   bmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nbin,npt,i,bin
      real pdata(npt),bdatax(nbin),bdatay(nbin),datamin,datamax,
     .   binsize,bmax,tsum

C     get size of bin.
      binsize=(datamax-datamin)/real(nbin-1)
      
C     initalize bdata
      do 20 i=1,nbin
         bdatay(i)=0.
 20   continue

C     loop over data and place data into proper bin.
      do 10 i=1,npt
C        calculate bin number
         bin=int((pdata(i)-datamin)/binsize)+1
         if((bin.gt.0).and.(bin.le.nbin)) bdatay(bin)=bdatay(bin)+1.0
 10   continue

C     get the max value in a bin and assign bin value.
      tsum=0.
      bmax=0.
      do 30 i=1,nbin
         bmax=max(bdatay(i),bmax)
         bdatax(i)=datamin+real(i-1)*binsize !added -1
         tsum=tsum+bdatay(i)
 30   continue
c      write(6,*) "Sum:",tsum
 
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m
      REAL*8 ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL*8 s,ep
      ave=0.0
      m=min(1000,n)
      do 11 j=1,m
         ave=ave+data(j)
 11   continue
      ave=ave/m
      var=0.0
      ep=0.0
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         var=var+s*s
 12   continue
      var=(var-ep**2/n)/(n-1) !Corrected two-pass formula (14.1.8).
      return
      END 