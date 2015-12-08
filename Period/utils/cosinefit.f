ccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cosinefit(pers,res,id,xcoo,ycoo,avg,a,btheta,rmavg,
     .     plot,ma,aerr,panx,pany,w,nstar,fixfreq)
CCCCccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     The one and only Fourier Decomposition routine.
C     per - fixed period for data
C     id,xcoo,ycoo - id and co-ordinates of object.
C     data is passed by common block /data/ 
C       x - time of observations
C       y - the observations
C       sig - the associated error
C       n - number of data points
      implicit none
      integer nmax,nfmax,nfit,nfitm,b,nstar,nstarmax,plot,endfit,panx,
     .  pany,fixfreq
      parameter(nmax=2000000,nfmax=2000,nfitm=2000,nstarmax=300)
      double precision w(nstarmax), x(nmax),y(nmax),sig(nmax),a(nfmax),
     .  chisq,tx,ty,ph(nmax),alpha(nfmax,nfmax),covar(nfmax,nfmax),
     .  alamda,avg,
     .     avgc,per,ty2(nmax),tx2(nmax),avgbin(10),maxy,miny,amp,bad,
     .     xcoo,ycoo,xi(nfitm*2),ctrl(12),f,pi,maxamp,
     .     fit1a,temp(2),integrate,btheta,rmavg,pers(nstarmax),arg,
     .     twopi,res(nmax),ochi,alamhi,chitol,aerr(nfmax)
      integer(kind=1) now(3)
      integer i,n,ia(nfmax),ma,nca, j, k,avgn(10),id,seed,status,cc,j1,
     .     i2,niter
      character*50 filename,starname
      logical loop
      common /data/ x,y,sig,n,nfit
c      external fit1a

      alamhi=100000.
      chitol=10e-9

      maxamp=2.0
      pi=3.141592654
      twopi=2*pi
      do 4 i=1,nstar
         w(i)=(2.0*pi)/pers(i)
c         write(6,*) "frequency",w(i),pers(i)
 4    continue

      loop=.true.
      i=0
      avg=0.
      maxy=y(1)
      miny=y(1)
      do 34 i=2,n
         maxy=max(maxy,y(i))
         miny=min(miny,y(i))
c         ph(i)=x(i)/pers(1) - int(x(i)/pers(1))
 34   continue
      
      avg=integrate(n,x,y,pers(1),0,1)      
      
      do 5 i=1,n
         y(i)=y(i)-avg
 5    continue
      
      ma=nstar*nfit*2+nstar
      do 10 i=1,ma
         ia(i)=1.0
 10   continue
      nca=nfmax

      alamda=-1.0

      if(xcoo.le.0.0)then
         j=(ma-nstar)/nstar
         do 30 k=nstar,nstar
            cc=(j+1)*(k-1)+2
            a(cc-1)=w(k)
            per=pers(k)
            do 31 i=cc,cc+j-2,2
               if(nstar.gt.1) then
                  temp(1)=integrate(n,x,res,per,(i-cc+2)/2,1)
                  temp(2)=integrate(n,x,res,per,(i-cc+2)/2,2)
               else
                  temp(1)=integrate(n,x,y,per,(i-cc+2)/2,1)
                  temp(2)=integrate(n,x,y,per,(i-cc+2)/2,2)
               endif
               a(i+1)=atan(-temp(1)/temp(2))
               a(i)=temp(2)/cos(a(i+1))
               a(i+1)=a(i+1)+pi
 31         continue
 30      continue

         do 32 k=1,nstar
            cc=(j+1)*(k-1)+2
C     THIS LINE HANDLES FIXED PERIODS
            if(fixfreq.gt.0) ia(cc-1)=0
            do 33 i=cc,cc+j-2,2
               if(a(i).lt.0.0) then
                  a(i)=abs(a(i))
                  a(i+1)=a(i+1)+pi
               endif
               if(a(i+1).lt.0) a(i+1)=2.0*pi+a(i+1)
               temp(1)=a(i+1)/(2.0*pi)
               b=int(temp(1))
               a(i+1)=a(i+1)-dble(b)*2.0*pi
 33         continue
 32      continue

      endif


c      goto 50

      endfit=0.
      i2=0
      niter=275
c      do 20 i2=1,niter
      do while (loop)
         i2=i2+1
         ochi=chisq
         call mrqmin(x,y,sig,n,a,ia,ma,covar,alpha,nca,chisq,alamda,w,
     .        nstar)
c         if(i2.ge.niter) then
         if((i2.gt.niter).or.(endfit.gt.50.0).or.(alamda.gt.alamhi))then
c            write(6,*) "i2,endfit,alamda,niter,alamhi"
c            write(6,*) "converge?",i2,endfit,alamda,niter,alamhi
            alamda=0.0
            call mrqmin(x,y,sig,n,a,ia,ma,covar,alpha,nca,chisq,alamda,
     .           w,nstar)
            loop=.false.
         endif

C     correct average
         avgc=0.
         do 35 j1=1,n
            ty=0.
            j=(ma-nstar)/nstar
            do 36  k=1,nstar
               cc=(j+1)*(k-1)+2
               do 37 i=cc,cc+j-2,2
                  arg=dble((i-cc+2)/2)*a(cc-1)*x(j1)+a(i+1)
                  ty=ty+a(i)*cos(arg)
 37            continue
 36         continue
            avgc=avgc+y(j1)-ty
 35      continue
         avgc=avgc/n
         do 16 j=1,n
            y(j)=y(j)-avgc
 16      continue
         avg=avg+avgc
C         write(6,*) avgc

C     look for negitive amplitudes!
CC         j=(ma-nstar)/nstar
CC         do 38 k=1,nstar
CC            cc=(j+1)*(k-1)+2
CC            do 39 i=cc,cc+j-2,2
CC               if(a(i).lt.0.0) then
CC                  a(i)=abs(a(i))
CC                  a(i+1)=a(i+1)+pi
CC               endif
CC               if (a(i+1).lt.0) a(i+1)=twopi+a(j+1)
CC               b = int(a(i+1)/twopi)
CC               a(i+1)=a(i+1)-b*twopi
CC 39         continue
CC 38      continue
         if(chisq.gt.ochi) then
            endfit=0
         elseif(abs(chisq-ochi).lt.chitol) then
            endfit=endfit+1
         endif
      enddo
c 20   continue
 50   continue
      niter=i2

      amp=0.

      do 21 i=1,n
         y(i)=y(i)+avg
 21   continue

      j=(ma-nstar)/nstar
      do 23 k=1,nstar
         cc=(j+1)*(k-1)+2
         a(cc-1)=abs(a(cc-1))
         pers(k)=twopi/a(cc-1)
         w(k)=1.0/pers(k)
 23   continue

      maxy=-99.9d10
      miny= 99.9d10

      do 22 i2=1,1000
         j=(ma-nstar)/nstar
         do 40 k=nstar,nstar
            cc=(j+1)*(k-1)+2
            per=twopi/a(cc-1)
            tx=dble(i2-1)*per/1000.0+x(1)
            tx2(i2)=tx/per-int(tx/per)
            ty2(i2)=0.
            do 41 i=cc,cc+j-2,2
               arg=real((i-cc+2)/2)*a(cc-1)*tx+a(i+1)
               ty2(i2)=ty2(i2)+a(i)*cos(arg)
 41         continue
 40      continue
         ty2(i2)=ty2(i2)+avg
         maxy=max(maxy,ty2(i2))
         miny=min(miny,ty2(i2))
 22   continue

      if(plot.eq.1) then
         if(nstar.eq.1) then 
            call plotph(n,x,y,per,panx,pany)
         else
            call plotph(n,x,res,per,panx,pany)
         endif
         call pgsci(2)
         call pgpt(1000,tx2,ty2,-1)
         call pgsci(1)
      endif
      amp = maxy-miny
c      if(1.0/per.gt.5.0) then
c      write(6,*) "freq,btheta,avg,amp,rchi,niter,alamda"
c      write(6,502)1.0/per,btheta,avg+rmavg,amp,chisq/(n-1),niter,alamda
 502  format(F9.4,1X,F10.3,1X,F9.4,1X,F9.4,1X,F9.4,1X,I4,F9.4)
      
      j=(ma-nstar)/nstar
      do 42 k=1,nstar
         cc=(j+1)*(k-1)+2
         aerr(cc-1)=sqrt(covar(cc-1,cc-1))/twopi
         if(cc+j-2.eq.nfmax)write(6,*) "Increase nfmax to: ",cc+j-2
c         write(6,501) "Fre", a(cc-1)/twopi,aerr(cc-1)
         do 43 i=cc,cc+j-2,2
            aerr(i)=sqrt(covar(i,i))
            aerr(i+1)=sqrt(covar(i+1,i+1))
 43      continue
c         write(6,500) "Amp", (a(i)  ,i=cc,cc+j-2,2)
c         write(6,500) "err", (sqrt(abs(covar(i,i))),    i=cc,cc+j-2,2)
c         write(6,500) "Psi", (a(i+1),i=cc,cc+j-2,2)
c         write(6,500) "err", (sqrt(abs(covar(i+1,i+1))),i=cc,cc+j-2,2)
 42   continue

c      endif
 500  format(A3,1X,20(F11.6))
 501  format(A3,1X,F11.5,F11.5)
C      CALL PGCLOS()
      rmavg=avg

      GOTO 999

C      close(10)
C      goto 999

 901  write(6,*) 'could not open ', filename
      goto 999
 902  write(6,*) 'error in file called' , filename
      goto 999

 999  return
      end
