CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program rhorand
C     generate synthetic rhoboot file
      implicit none
      integer now(3),seed,i,niter,np,j,accrho
      parameter(np=7)
      real rhoi,rhoierr(9),gasdev,ran2,dumr,out(np),rhoin(9),rgas,rnorm
      real M,Mpe,Mme,R,Rpe,Rme,Rsun,Msun,Pi,drho,dsig,ranrho
c      data rhoin/-1.0d2,-3.0d0,-2.0d0,-1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,
c     .  1.0d2/
      data rhoin/-0.99999999,-0.9973,-0.954,-0.683,0.0d0,0.683,0.954
     .   ,0.9973,0.99999999/


      Msun=1.989e30
      Rsun=695500000.0
      Pi=acos(-1.d0)!define Pi

C     initialize random number generator
      call itime(now)
      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      dumr=ran2(-seed)

!     rhostar from transit models.
      rhoi=8.6
      rhoierr(1)=-rhoi
      rhoierr(2)=-3.3
      rhoierr(3)=-2.2
      rhoierr(4)=-1.1
      rhoierr(5)=0.0d0
      rhoierr(6)= 1.2
      rhoierr(7)= 2.4
      rhoierr(8)= 3.6
      rhoierr(9)=rhoierr(8)*10.0d0

!     Mass+rad
      M=0.51
      Mpe= 0.06
      Mme=-0.06
      R=0.4476
      Rpe= 0.0275
      Rme=-0.0275

!     rhostar from transit models.
      rhoi=0.837
      rhoierr(1)=-rhoi
      rhoierr(2)=-0.564
      rhoierr(3)=-0.376
      rhoierr(4)=-0.188
      rhoierr(5)=0.0d0
      rhoierr(6)= 0.396
      rhoierr(7)= 0.792
      rhoierr(8)= 1.188
      rhoierr(9)=rhoierr(8)*10.0d0

!     Mass+rad
      M=1.037
      Mpe= 0.054
      Mme=-0.047
      R=1.109
      Rpe= 0.147
      Rme=-0.091

!     rhostar from transit models.
      rhoi=0.837
      rhoierr(1)=-rhoi
      rhoierr(2)=-0.564
      rhoierr(3)=-0.376
      rhoierr(4)=-0.188
      rhoierr(5)=0.0d0
      rhoierr(6)= 0.396
      rhoierr(7)= 0.792
      rhoierr(8)= 1.188
      rhoierr(9)=rhoierr(8)*10.0d0

!     Mass+rad
      M=0.20
      Mpe= 0.04
      Mme=-0.04
      R=0.247
      Rpe= 0.032
      Rme=-0.032

      niter=100000
      do 12 j=1,niter

         rgas=gasdev(seed)
         out(2)=6.0+2.0*rgas
c         out(2)=7.97+2.0*rgas
         rgas=gasdev(seed)
c         out(3)=-0.28+0.10*rgas
c         out(3)=0.21+0.09*rgas
         out(3)=-0.2+0.11*rgas
         rgas=gasdev(seed)
c         out(6)=3841.0+50.0*rgas
c         out(6)=5757.0+85.0*rgas
         out(6)=3244.0+44.0*rgas

         accrho=0
         do while(accrho.eq.0)
            rgas=gasdev(seed)
            if(rgas.ge.0.0d0)then
               out(1)=M+Mpe*rgas
            else
               out(1)=M-Mme*rgas
            endif

            rnorm=ran2(seed)
            rgas=gasdev(seed)
c            if(rnorm.lt.abs(Rme/Rpe))then
c               rgas=abs(rgas)
c            else
c               rgas=-abs(rgas)
c            endif
            if(rgas.ge.0.0d0)then
               out(4)=R+Rpe*rgas
            else
               out(4)=R-Rme*rgas
            endif

            out(5)=out(1)*Msun/(4.0/3.0*pi*(Rsun*out(4))**3.0)

c            drho=out(5)/1000.0d0-rhoi
c            call getrhosig(rhoierr,rhoin,9,drho,dsig)
c            ranrho=ran2(seed)
c            if(ranrho.gt.abs(dsig))then
               accrho=1
c            else
c               accrho=0
c            endif

!            write(0,*) ranrho,abs(dsig),out(5)/1000.0d0
!            write(0,*) "accrho: ",accrho
!            read(5,*)
         enddo

         out(7)=log10(out(4)**2.0*(out(6)/5781.6)**4.0)

         write(6,501) (out(i),i=1,np),0.0,0,1
 501     format(8(1X,1PE17.10),2(1X,I2))
 12   continue

 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getrhosig(rhoierr,rhoin,npt,drho,dsig)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      real rhoierr(npt),rhoin(npt),drho,dsig

C     Default is to reject drho is not within bounds below
      dsig=100.0d0

c      write(0,*) (rhoin(i),i=1,npt)
c      write(0,*) (rhoierr(i),i=1,npt)

c      write(0,*) drho,rhoierr(2),rhoierr(8)
      if(drho.lt.rhoierr(2))then
        dsig=rhoin(2)+(drho-rhoierr(2))/(rhoierr(3)-rhoierr(2))*
     .          (rhoin(3)-rhoin(2))
      elseif(drho.gt.rhoierr(8))then
        dsig=rhoin(8)+(drho-rhoierr(8))/(rhoierr(9)-rhoierr(8))*
     .          (rhoin(9)-rhoin(8))
      else
        do 10 i=2,npt-2
            if((drho.gt.rhoierr(i)).and.(drho.le.rhoierr(i+1)))then
                dsig=rhoin(i)+(drho-rhoierr(i))/(rhoierr(i+1)-
     .              rhoierr(i))*(rhoin(i+1)-rhoin(i))
            endif
 10     continue
      endif
c      write(0,*) "dhro",drho,dsig
c      read(5,*)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getrhosigold(rhoierr,rhoin,npt,drho,dsig)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer npt,i
      real rhoierr(npt),rhoin(npt),drho,dsig

C     Default is to reject drho is not within bounds below
      dsig=100.0d0

c      write(0,*) (rhoin(i),i=1,npt)
c      write(0,*) (rhoierr(i),i=1,npt)
      do 10 i=1,npt-1
        if((drho.gt.rhoierr(i)).and.(drho.le.rhoierr(i+1)))then
            dsig=rhoin(i)+(drho-rhoierr(i))/(rhoierr(i+1)-rhoierr(i))*
     .          (rhoin(i+1)-rhoin(i))
        endif
 10   continue
c      write(0,*) "dhro",drho,dsig
c      read(5,*)

      return
      end


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
 1       v1=2.*ran2(idum)-1.
         v2=2.*ran2(idum)-1.
         rsq=v1**2+v2**2
         if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
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
      real FUNCTION ran2(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL AM,EPS,RNMX
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
