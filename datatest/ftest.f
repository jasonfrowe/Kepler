      program ftest
      implicit none
      integer nmax,n1,n2,i,iargc,nunit
      parameter(nmax=600000)
      double precision data1(nmax),data2(nmax),f,prob,dumr
      character*80 file1,file2
      
      if(iargc().lt.2) goto 903
      
      call getarg(1,file1)
      call getarg(2,file2)
      
      nunit=10
      open(unit=nunit,file=file1,status='old',err=901)
      i=1
 10   read(nunit,*,end=11) dumr,data1(i)
        write(0,*) i,data1(i)
        i=i+1
      goto 10
 11   continue
      n1=i-1
      close(nunit)
      
      open(unit=nunit,file=file2,status='old',err=902)
      i=1
 12   read(nunit,*,end=13) dumr,data2(i)
        write(0,*) i,data2(i)
        i=i+1
      goto 12
 13   continue
      n2=i-1
      close(nunit)
      write(0,*) "n1,n2:",n1,n2
      
      call ftestnr(data1,n1,data2,n2,f,prob)
      write(0,*) "F,Prob:",f,prob
      
      goto 999
 901  write(0,*) "Cannot open ",file1
      goto 999
 902  write(0,*) "Cannot open ",file2
      goto 999
 903  write(0,*) "Usage: file1 file2"
      write(0,*) "column 2 of fileN is used to calculate variances"
      write(0,*) "Output: f prob"
      write(0,*) "Small values of prob indicate that the two arrays have
     . significantly different variances."
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ftestnr(data1,n1,data2,n2,f,prob)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n1,n2
      REAL*8 f,prob,data1(n1),data2(n2)
C     USES avevar,betai
C     Given the arrays data1(1:n1) and data2(1:n2), this routine returns 
C     the value of f, and its significance as prob. Small values of prob
C     indicate that the two arrays have significantly different 
C     variances.
      REAL*8 ave1,ave2,df1,df2,var1,var2,betai
      call avevar(data1,n1,ave1,var1)
      call avevar(data2,n2,ave2,var2)
C     Make F the ratio of the larger variance to the smaller one.
      if(var1.gt.var2)then 
        f=var1/var2
        df1=n1-1
        df2=n2-1
      else
        f=var2/var1
        df1=n2-1
        df2=n1-1
      endif
      prob=2.*betai(0.5*df2,0.5*df1,df2/(df2+df1*f))
      if(prob.gt.1.)prob=2.-prob
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n
      REAL*8 ave,var,data(n)
C     Given array data(1:n), returns its mean as ave and its variance 
C     as var.
      INTEGER j
      REAL*8 s,ep
      ave=0.0
      do 11 j=1,n
        ave=ave+data(j)
 11   continue
      ave=ave/n
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision FUNCTION betai(a,b,x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      REAL*8 a,b,x
C     USES betacf,gammln
C     Returns the incomplete beta function Ix(a, b).
      REAL*8 bt,betacf,gammln
      if(x.lt.0..or.x.gt.1.)pause "bad argument x in betai"
      if(x.eq.0..or.x.eq.1.)then
        bt=0.
      else !Factors in front of the continued fraction.
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)
     *      +a*log(x)+b*log(1.-x))
      endif
      if(x.lt.(a+1.)/(a+b+2.))then !Use continued fraction directly.
        betai=bt*betacf(a,b,x)/a
        return
      else
C     Use continued fraction after making the symmetry transformation.
        betai=1.-bt*betacf(b,a,1.-x)/b 
        return
      endif
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision FUNCTION betacf(a,b,x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER MAXIT
      REAL*8 a,b,x,EPS,FPMIN
      PARAMETER (MAXIT=100,EPS=3.e-7,FPMIN=1.e-30)
C     Used by betai: Evaluates continued fraction for incomplete beta 
C     function by modified Lentz’s method (§5.2).
      INTEGER m,m2
      REAL*8 aa,c,d,del,h,qab,qam,qap
      qab=a+b  !These q’s will be used in factors that occur in the
      qap=a+1. !coefficients (6.4.6).
      qam=a-1.
      c=1.     !First step of Lentz’s method.
      d=1.-qab*x/qap
      if(abs(d).lt.FPMIN)d=FPMIN
      d=1./d
      h=d
      do 11 m=1,MAXIT
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.+aa*d !One step (the even one) of the recurrence.
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        h=h*d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.+aa*d !Next step of the recurrence (the odd one).
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1 !Are we done?
 11   continue
      pause "a or b too big, or MAXIT too small in betacf"
 1    betacf=h
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION gammln(xx)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      REAL*8 xx
C     Returns the value ln[Γ(xx)] for xx > 0.
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
C     Internal arithmetic will be done in double precision, a nicety 
C     that you can omit if five-figure accuracy is good enough.
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     * 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     * -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
 11   continue
      gammln=tmp+log(stp*ser/x)
      return
      END