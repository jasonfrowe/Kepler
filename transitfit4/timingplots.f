      program timingplot
      implicit none
      integer i,nmax,npt,nunit,nunit2,j,npt2,now(3),seed,nfp,kid,k
      parameter(nmax=10000)
      integer fp(nmax),ns(nmax)
      real x(nmax),y(nmax)
      double precision bjd(nmax),oc(nmax),ocerr(nmax),koinum,chisq,
     .  chi(nmax),rchi(nmax),dumr,ran2,gasdev,sigma(nmax),koi(nmax),
     .  sn(nmax),snr,medocerr,absoc(nmax),medabsoc
      character*80 filename,filelist,fplist

      call itime(now)
      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      dumr=ran2(-seed)

      nunit=10
      fplist="koi_FPs.20120229.csv"
      i=1
      open(unit=nunit,file=fplist,status='old',err=903)
 5    read(nunit,*,end=6) koi(i),kid,sn(i),fp(i)
        if(koi(i).ge.2669.0) fp(i)=1 !ignore Q1-Q8 KOIs
        if(koi(i)-int(koi(i)).eq.0.0) koi(i)=koi(i)+0.01
c        if(fp(i).eq.1)then
cc            write(0,*) koi(i),kid,fp(i)
c            i=i+1
c        endif
        i=i+1
      goto 5
 6    continue
      nfp=i-1
      write(0,*) "nFP: ",nfp
      close(nunit)

      
      filelist="ttplots.list"
      nunit=10
      nunit2=11
      j=1
      open(unit=nunit,file=filelist,status='old',err=901)
 10   read(nunit,*,end=11) filename
        read(filename(4:10),*) koinum
        
        do 16 i=1,nfp
            if(abs(koinum-koi(i)).lt.0.005)then
c                if(fp(i).eq.1)goto 10
                snr=sn(i)
            endif
 16     continue
 
        open(unit=nunit2,file=filename,status='old',err=902)

        i=1
 12     read(nunit2,*,end=13) bjd(i),oc(i),ocerr(i)
            if(ocerr(i).eq.0.0d0) goto 12
            ocerr(i)=ocerr(i)*1.18d0
            i=i+1
        goto 12
 13     continue
        npt=i-1
        
C       Remove points with ocerr > 2 * median(ocerr)
        call rqsort(npt,ocerr,ns)
        medocerr=ocerr(ns(npt/2+1))       
        k=0
        do 17 i=1,npt
            if(ocerr(i).lt.2.0d0*medocerr)then
                k=k+1
            endif
            bjd(k)=bjd(i)
            oc(k)=oc(i)
            ocerr(k)=ocerr(i)
 17     continue
        npt=k

C       Remove points with abs(TTV) > 4 x median(abs(TTV))
        do 18 i=1,npt
            absoc(i)=abs(oc(i))
 18     continue
        call rqsort(npt,absoc,ns)
        medabsoc=absoc(ns(npt/2+1))
        k=0
        do 19 i=1,npt
            if(absoc(i).lt.4.0d0*medabsoc)then
                k=k+1
            endif
            bjd(k)=bjd(i)
            oc(k)=oc(i)
            ocerr(k)=ocerr(i)
 19     continue
        npt=k
        
        if(npt.lt.1) goto 10 !make sure we have at least one point
        
        snr=snr/sqrt(dble(npt))
        chi(j)=chisq(npt,oc,ocerr)
        rchi(j)=chi(j)/dble(npt-1)
        sigma(j)=sqrt(2.0d0/dble(npt-1))
 
        write(6,500) koinum,rchi(j),(rchi(j)-1.0)/sigma(j),snr,npt
c        read(5,*)
 500    format(F7.2,' &',3(1X,F7.2,' &'),1X,I4,' \\')      
        
        close(nunit2)

        j=j+1
      goto 10
 11   continue     
      close(nunit)
      npt2=j-1
      
      do 14 i=1,npt2
        x(i)=real(rchi(i))
 14   continue

      call pgopen('?')
      call PGPAP ( 7.0 ,1.0)  
            
      do 15 i=1,npt2
        y(i)=real(1.0+sigma(i)*gasdev(seed))
 15   continue
 
      call PGHIST(npt2,y,0.0,20.0,50,0)
      call pgsci(2)
      call PGHIST(npt2,x,0.0,20.0,50,1)
      call pgsci(1)
      call pglabel("reduced chi-squared","Frequency (N)","")
    
      call PGHIST(npt2,y,0.0,5.0,50,0)
      call pgsci(2)
      call PGHIST(npt2,x,0.0,5.0,50,1)
      call pgsci(1)  
      call pglabel("reduced chi-squared","Frequency (N)","")

      
      call pgclos()
      
      goto 999
 901  write(0,*) "Cannot open ",filelist
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 903  write(0,*) "Cannot open ",fplist
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function chisq(npt,oc,ocerr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      double precision oc(npt),ocerr(npt)
      
      chisq=0.0d0
      do 10 i=1,npt
        chisq=chisq+oc(i)*oc(i)/(ocerr(i)*ocerr(i))
 10   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision FUNCTION ran2(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     * IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     * IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.-EPS)
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
      double precision FUNCTION gasdev(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,ran2
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
      
c**********************************************************************
      subroutine rqsort(n,a,p)
c======================================================================
c     Return integer array p which indexes array a in increasing order.
c     Array a is not disturbed.  The Quicksort algorithm is used.
c
c     B. G. Knapp, 86/12/23
c
c     Reference: N. Wirth, Algorithms and Data Structures,
c     Prentice-Hall, 1986
c======================================================================
      implicit none

c     Input:
      integer   n
      real*8      a(n)

c     Output:
      integer   p(n)

c     Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
c        (LGN = log base 2 of maximum n;
c         Q = smallest subfile to use quicksort on)

c     Local:
      real*8      x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

c     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

c     Initialize the pointer array
      do 1 i=1,n
         p(i)=i
    1 continue

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

c           Use straight insertion
            do 6 i=l+1,r
               t = p(i)
               x = a(t)
               do 4 j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
    4          continue
               j=l-1
    5          p(j+1) = t
    6       continue
         else

c           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

c           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

c           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end