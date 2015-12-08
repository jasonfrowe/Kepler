      subroutine polyfilter(npt,time,flux,ferr,ts,tflag,boxbin,x,y,z,
     .  ngap,gaps,offset,nfitp,work)
      implicit none
      integer npt,ts(npt),tflag(npt),i,j,k,ngap,nfitp,npt2
      double precision time(npt),flux(npt),ferr(npt),tzero,boxbin,bbd2,
     .  x(npt),y(npt),z(npt),t1,t2,gaps(npt),cadence,nc,offset(npt),off,
     .  mean,mx,work(npt)

      bbd2=boxbin/2.0d0

C     remove 30 for LC data.
      nc=3000.0d0 !number of cadences that need to be missing to mark a gap
      cadence=0.02*nc

C     sort the data by time
c      write(0,*) "sorting"
      call rqsort(npt,time,ts)
c      write(0,*) "done sorting"
      
      ngap=0
      do 12 i=1,npt-1
        if(time(ts(i+1))-time(ts(i)).gt.cadence)then
            ngap=ngap+1
            gaps(ngap)=(time(ts(i+1))+time(ts(i+1)))/2.0d0
c            write(0,*) ngap,gaps(ngap)
        endif
 12   continue
        
      
      
      do 10 i=1,npt
c      do 10 i=48600,npt
c        write(0,*) "i:",i
        tzero=time(ts(i))
        t1=tzero-bbd2
        t2=tzero+bbd2
        do 13 j=1,ngap
            if((gaps(j).gt.t1).and.(gaps(j).lt.tzero))t1=gaps(j)
            if((gaps(j).lt.t2).and.(gaps(j).gt.tzero))t2=gaps(j)
 13     continue
        k=0 !reset counter
        do 11 j=1,npt
            if((time(j).gt.t1).and.(time(j).lt.t2).and.
     .        (tflag(j).eq.0))then
                k=k+1
                x(k)=time(j)
                y(k)=flux(j)
                z(k)=ferr(j)
c                write(6,*) x(k),y(k),z(k)
            endif
 11     continue
        npt2=k
c        std=stdev(npt2,x,mx)
c        do 
c        write(0,*) t1,t2,k
        if(npt2.gt.nfitp+1)then
c            write(0,*) "into polydetrend"
            call polydetrend(npt2,x,y,z,nfitp,tzero,off,work)
c            write(0,*) "out of polydetrend"
            offset(i)=off
        elseif(npt2.eq.0)then
            offset(i)=0.0d0
        else
            offset(i)=mean(npt2,y)
        endif
c        write(0,*) "off",off
c        read(5,*)
 10   continue
      
      do 14 i=1,npt
        flux(ts(i))=flux(ts(i))-offset(i)
 14   continue
        
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function mean(npt,pts)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      double precision pts(npt)
      
      mean=0.
      do 10 i=1,npt
        mean=mean+pts(i)
 10   continue
      mean=mean/dble(npt)
      if(npt.le.0) mean=0.
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine polydetrend(npt,time,mag,merr,nfit,tzero,off,work)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nfit
      integer i,ma,j,npc,ii,nmax,maxiter
      parameter(npc=10,nmax=1800000,maxiter=10)
      integer ia(npc)
      real*8 time(npt),mag(npt),merr(npt),ans(npc),covar(npc,npc),
     .     chisq,T,meanT,off,tzero,work(npt),std,mx,stdev,sigcut,
     .   merr2(nmax),dchi,ochi
     
      sigcut=3.0

      if(nfit.gt.npc) pause "increase npc in detrend"
      do 36 i=1,nfit
         ia(i)=1
 36   continue
 
      meanT=0.0d0
      do 11 i=1,npt
        meanT=meanT+time(i)
 11   continue
      meanT=meanT/dble(npt)

      do 12 i=1,npt
        time(i)=time(i)-meanT
 12   continue

      ma=nfit
      call lfit(time,mag,merr,npt,ans,ia,ma,covar,npc,chisq)
c      write(0,*) (ans(i),i=1,nfit)

      ii=0
      dchi=1
      ochi=chisq
      do 14 while((ii.lt.maxiter).and.(dchi.gt.0.1))
c         write(0,*) "ii:",ii
!      do 14 ii=1,10
!         write(0,*) "npt: ",npt

         do 10 i=1,npt
            T=0.
            DO 33 J = 1, NFIT
               T = ANS(J)*((time(i))**REAL(J-1)) +  T
 33         CONTINUE
            work(i)=mag(i)-T
 10      continue

         std=stdev(npt,work,mx)
      
         j=0
         do 13 i=1,npt
            if(abs(work(i)-mx).lt.sigcut*std)then
               j=j+1
               time(j)=time(i)
               mag(j)=mag(i)
               merr2(j)=merr(i)+abs(work(i))
            endif
 13      continue
         npt=j
         call lfit(time,mag,merr2,npt,ans,ia,ma,covar,npc,chisq)
         dchi=abs(chisq-ochi)
         ochi=chisq
!         write(0,*) "npt: ",npt
         ii=ii+1  !count number of iterations to avoid infinite loops
 14   continue

      T=0.
      DO 34 J = 1, NFIT
         T = ANS(J)*((tzero-meanT)**REAL(J-1)) +  T
 34   CONTINUE
      off=T

c      do 13 i=1,npt
c        time(i)=time(i)+meanT
c 13   continue

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE lfit(x,y,sig,ndat,a,ia,ma,covar,npc,chisq) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER ma,ia(ma),npc,ndat,MMAX
      REAL*8 chisq,a(ma),covar(npc,npc),sig(ndat),x(ndat),y(ndat)
      EXTERNAL funcs2   
      PARAMETER (MMAX=50) !Set to the maximum number of coefficients ma.
C     USES covsrt,gaussj
      INTEGER i,j,k,l,m,mfit
      REAL*8 sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX)

      mfit=0
      do 11 j=1,ma
        if(ia(j).ne.0) mfit=mfit+1
 11   continue
      if(mfit.eq.0) pause "lfit: no parameters to be fitted"
      do 13 j=1,mfit !Initialize the (symmetric) matrix.
        do 12 k=1,mfit
            covar(j,k)=0.
 12     continue
        beta(j)=0.
 13   continue
      do 17 i=1,ndat !Loop over data to accumulate coefficients of the normal
        call funcs2(x(i),afunc,ma) !equations.
        ym=y(i)
        if(mfit.lt.ma) then !Subtract off dependences on known pieces of the fitting
            do 14 j=1,ma !function.
                if(ia(j).eq.0) ym=ym-a(j)*afunc(j)
 14         continue
        endif
        sig2i=1./sig(i)**2
        j=0
        do 16 l=1,ma
            if (ia(l).ne.0) then
                j=j+1
                wt=afunc(l)*sig2i
                k=0
                do 15 m=1,l
                    if (ia(m).ne.0) then
                        k=k+1
                        covar(j,k)=covar(j,k)+wt*afunc(m)
                    endif
 15             continue
                beta(j)=beta(j)+ym*wt
            endif
 16     continue
 17   continue
      do 19 j=2,mfit !Fill in above the diagonal from symmetry.
        do 18 k=1,j-1
            covar(k,j)=covar(j,k)
 18     continue
 19   continue
      call gaussj(covar,mfit,npc,beta,1,1) !Matrix solution.
      j=0
      do 21 l=1,ma
        if(ia(l).ne.0) then
            j=j+1
            a(l)=beta(j) !Partition solution to appropriate coefficients a.
        endif
 21   continue
      chisq=0.
      do 23 i=1,ndat
        call funcs2(x(i),afunc,ma)
        sum=0.
        do 22 j=1,ma
            sum=sum+a(j)*afunc(j)
 22     continue
        chisq=chisq+((y(i)-sum)/sig(i))**2
 23   continue
      call covsrt(covar,npc,ma,ia,mfit) !Sort covariance matrix to true order of fitting
      return
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FUNCS2(X,P,NP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
C
C     Taken from NUMERICAL RECIPES: The Art of Scientific Computing by
C     Press, Flannery, Teukolsky, Vetterling (1986).
C
      integer j,np
      real*8 P(NP),X
      P(1)=1.
      DO 11 J=2,NP
         P(J)=P(J-1)*X
 11   CONTINUE
      RETURN
      END      
