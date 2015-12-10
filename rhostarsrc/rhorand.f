CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program rhorand
C     generate synthetic rhoboot file
      implicit none
      integer iargc,now(3),seed,i,niter,np,j
      parameter(np=14)
      real gasdev,gas(np),out(np/2),dumr,ran2
      character*80 cline

C     initialize random number generator
      call itime(now)
      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      dumr=ran2(-seed)
      
      
      if(iargc().lt.15) goto 901

      do 10 i=1,np
        call getarg(i,cline)
        read(cline,*) gas(i)
 10   continue
      call getarg(15,cline)
      read(cline,*) niter

      do 12 j=1,niter
         do 11 i=1,np,2
c         write(0,*) i,i+1,i/2+1
            out(i/2+1)=gasdev(seed)*gas(i+1)+gas(i)
 11      continue
         out(5)=out(5)*1000.0 !convert g/cm^3 to kg/m^3
         if(out(7).lt.0.0) out(7)=0.1
         out(7)=log10(out(7))
         write(6,501) (out(i),i=1,np/2),0.0,0, 1
 501     format(8(1X,1PE17.10),2(1X,I2))
 12   continue
 
      goto 999
 901  write(0,*) "Usage: M Merr Age Ageerr Z Zerr R Rerr rhostar rhostar
     .err Teff Tefferr L Lerr niter"
      write(0,*) " niter: number of points to generate"
      goto 999
 999  end
      
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
