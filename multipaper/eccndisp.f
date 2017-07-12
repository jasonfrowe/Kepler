      program eccndisp
      implicit none
      integer nmax,nfit,i,j,flag,nc,nmc,nunit,npt,neccn,now(3),seed
      parameter(nmax=500000,nfit=5,neccn=100)
      integer nbstar(nmax)
      real chi(nmax),ans(nfit,nmax),Pbin(neccn,nmax),Psum(neccn),deccn,
     .   eccn,Pbeta,a,b,eccnp(neccn),bx(nmax),ave,var,Perr(neccn),dumr,
     .   ran2,gasdev,Psum2(neccn),Perr2(neccn),chisq,xp(2*neccn),
     .   yp(2*neccn)
      character*80 mcmcfile

      if(iargc().lt.1) goto 901

      call getarg(1,mcmcfile)

      nunit=10
      open(unit=nunit,file=mcmcfile,status='old',err=902)
      i=1
 10   read(nunit,500,end=11) flag,nc,nmc,chi(i),(ans(j,i),j=1,nfit),
     . nbstar(i)
 500     format(I1,2(1X,I1),6(1X,1PE17.10),1X,I5)
         i=i+1
      goto 10
 11   continue
      close(nunit)
      npt=i-1
      write(0,*) "npt: ",npt

      do 15 i=1,npt
         do 14 j=1,neccn
            Pbin(j,i)=0.0
 14      continue
 15   continue

      deccn=1.0/real(neccn+1)
      do 12 i=1,npt
         do 13 j=1,neccn
            eccn=deccn*real(j)
            a=ans(2,i)
            b=ans(3,i)
            Pbin(j,i)=Pbeta(eccn,a,b)
 13      continue
 12   continue

      do 16 j=1,neccn
         Psum(j)=0.0
 16   continue
      do 17 j=1,neccn
         do 18 i=1,npt
            bx(i)=Pbin(j,i)
c            Psum(j)=Psum(j)+Pbin(j,i)
 18      continue
         call avevar(bx,npt,ave,var)
         Psum(j)=ave
         Perr(j)=sqrt(var)
         eccnp(j)=deccn*real(j)
         xp(j)=eccnp(j)
         xp(2*neccn-j+1)=eccnp(j)
         yp(j)=ave+perr(j)
         yp(2*neccn-j+1)=ave-perr(j)
c         write(0,*) eccnp(j),Psum(j),Perr(j)
 17   continue
c      do 19 j=1,neccn
c         Psum(j)=Psum(j)/real(npt)
c         eccnp(j)=deccn*real(j)
c         write(6,*) eccnp(j),Psum(j)
c 19   continue

      call pgopen('?')
      call pgpage()
      call PGPAP ( 8.0 ,1.0)
      call pgsch(1.5)
      call pgslw(2)

      call pgvport(0.15,0.85,0.15,0.85)
      call pgwindow(0.0,1.0,0.0,4.0)
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Eccentricity","Probability Density","")
c      call pgline(neccn,eccnp,Psum)
c      call pgpt(neccn,eccnp,Psum,17)
c      call pgerrb(6,neccn,eccnp,Psum,Perr,1.0)
      call pgsci(15)
      call pgpoly(2*neccn,xp,yp)
      call pgsci(1)
      call pgline(neccn,eccnp,Psum,17)
      do i=1,neccn
         write(6,*) eccnp(i),Psum(i)
      enddo

      call itime(now)
      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      dumr=ran2(-seed)

      do 21 j=1,neccn
         do 20 i=1,npt
            eccn=deccn*real(j)
            a=0.867+0.044*gasdev(seed)
            b=3.030+0.165*gasdev(seed)
            bx(i)=Pbeta(eccn,a,b)
 20      continue
         call avevar(bx,npt,ave,var)
         Psum2(j)=ave
         Perr2(j)=sqrt(var)
         eccnp(j)=deccn*real(j)
         xp(j)=eccnp(j)
         xp(2*neccn-j+1)=eccnp(j)
         yp(j)=ave+perr2(j)
         yp(2*neccn-j+1)=ave-perr2(j)
c         write(0,*) eccnp(j),Psum2(j),Perr2(j)
 21   continue

c      call pgsci(2)
c     call pgpt(neccn,eccnp,Psum2,17)
c      call pgerrb(6,neccn,eccnp,Psum2,Perr2,1.0)
      call pgsci(5)
      call pgpoly(2*neccn,xp,yp)
      call pgsci(4)
      call pgline(neccn,eccnp,Psum2,17)


      call pgsci(1)

      call pgclos()

      chisq=0.0
      do 22 i=1,neccn
         chisq=chisq+(Psum(i)-Psum2(i))**2.0/(Perr(i)*Perr(i)+
     .    Perr2(i)*Perr2(i))
 22   continue
      write(6,*) "chisq: ",chisq



      goto 999
 901  write(0,*) "Usage: tdurplot mcmcfile"
      goto 999
 902  write(0,*) "Cannot open ",mcmcfile
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION ran2(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     * IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     * IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION gasdev(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
 1       v1=2.0d0*ran2(idum)-1.0d0
         v2=2.0d0*ran2(idum)-1.0d0
         rsq=v1**2+v2**2
         if(rsq.ge.1.0d0.or.rsq.eq.0.0d0)goto 1
         fac=sqrt(-2.0d0*log(rsq)/rsq)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      endif
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m
      REAL ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL s,ep
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function Pbeta(x,a,b)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real x,a,b,gammln,beta

      beta=exp(gammln(a)+gammln(b)-gammln(a+b))

      Pbeta=(1.0/beta)*x**(a-1.0)*(1-x)**(b-1.0)
c      write(6,*) Pbeta,a,b,gammln(a+b)
c      read(5,*)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION gammln(xx)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
C     From Numerical Recipes
      real xx
C     Returns the value ln[Γ(xx)] for xx > 0.
      INTEGER j
      double precision ser,stp,tmp,x,y,cof(6)
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
