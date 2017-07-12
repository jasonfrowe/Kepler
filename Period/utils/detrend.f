CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine detrend(npt,time,mag,merr,nfit)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nfit
      integer i,ma,j,npc
      parameter(npc=10)
      integer ia(npc)
      real*8 time(npt),mag(npt),merr(npt),ans(npc),covar(npc,npc),
     .     chisq,T

      if(nfit.gt.npc) pause "increase npc in detrend"
      do 36 i=1,nfit
         ia(i)=1
 36   continue

      ma=nfit
      call lfit(time,mag,merr,npt,ans,ia,ma,covar,npc,chisq)
      write(6,*) (ans(i),i=1,nfit)

      do 10 i=1,npt
         T=0.
         DO 33 J = 1, NFIT
            T = ANS(J)*((time(i))**REAL(J-1)) +  T
 33      CONTINUE
         mag(i)=mag(i)-T
 10   continue

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