C     Estimate star noise limits from CoRoT data.
C     Needs the file 'bstars.dat'
      program bootnoise
      implicit none
      integer nmax,i,nr,ns,now(3),seed,nboot,nG,nK,nM,ndAb,ndAf,
     .	ndG,ndK,ndM
      parameter(nmax=65,nboot=1000)
      real rmag(nmax),noise(nmax),ref(nmax),samp(nmax),ran2,dumr,mnoise,
     .	gnoise,merr,boottemp(nboot),nsamp(nmax),meanboot,sigma,
     .  dataAb(nmax),dataAf(nmax),dataG(nmax),dataK(nmax),datam(nmax),
     .	d,prob
      character*2 sptype(nmax)
      character*80 filename

      sigma=1.0 !number of sigma for "quiet" cut
     
      filename="bstars.dat"
      open(unit=10,file=filename,status='old',err=901)

C     Read in the data file.      
      do 10 i=1,nmax
      	read(10,*) rmag(i),noise(i),sptype(i)
 10   continue
      
C     Start with comparision of bright A-stars to G-stars.      
      
      nr=0 !initialize counter to track sample sizes.
      ns=0 
      ndAb=0
      ndG=0
      do 11 i=1,nmax
      	if((sptype(i)(1:1).eq.'A').and.(rmag(i).lt.13.0))then
      		nr=nr+1
      		ndAb=ndAb+1
      		ref(nr)=noise(i) !store A star noises
      		dataAb(ndAb)=noise(i)
        endif
        if(sptype(i)(1:1).eq.'G')then
        	ns=ns+1
        	ndG=ndG+1
        	samp(ns)=noise(i) !store G star noises
        	dataG(ndG)=noise(i)
        endif
 11   continue
      write(0,*) "Number of A-stars: ",nr
      write(0,*) "Number of G-stars: ",ns

      call kstwo(dataAb,ndAb,dataG,ndG,d,prob)
      write(0,*) "KS_AG d:",d," prob:",prob

C     Calculate mean noise for bright A-stars
      mnoise=0.0
      do 12 i=1,nr
      	mnoise=mnoise+ref(i)
 12   continue
      mnoise=mnoise/real(nr)
c      write(0,*) "A-star bright noise: ",mnoise

C     Calculate mean noise for G-stars
      gnoise=0.0
      do 13 i=1,ns
      	gnoise=gnoise+samp(i)
 13   continue
      gnoise=gnoise/real(ns)
c      write(0,*) "G-star bright noise: ",gnoise

C     Random number creation (so we get a different seed each time)
      call itime(now)
      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      dumr=ran2(-seed)

C     get the error on the mean from a simple bootstrap    
c      merr=meanboot(seed,ns,samp,nboot,boottemp,nmax,nsamp)
	  call avevar(ref,nr,mnoise,merr)
      merr=sqrt(merr)
      write(0,*) "A-star bright noise: ",mnoise,"+/-",merr
      
      nG=0 !number of G-stars that are "quiet"
      do 14 i=1,ns
      	if(samp(i).lt.mnoise+merr*sigma)then
      		nG=nG+1
        endif
 14   continue
 
      write(0,*) "Percentage of Quiet G-stars:",100.0*real(nG)/real(ns) 
      	
      ns=0
      ndK=0
      do 15 i=1,nmax
        if(sptype(i)(1:1).eq.'K')then
        	ns=ns+1
        	ndK=ndK+1
        	samp(ns)=noise(i) !store K star noises
        	dataK(ndk)=noise(i)
        endif
 15   continue
      write(0,*) "Number of K-stars: ",ns     
      
      call kstwo(dataAb,ndAb,dataK,ndK,d,prob)
      write(0,*) "KS_AK d:",d," prob:",prob
      call kstwo(dataG,ndG,dataK,ndK,d,prob)
      write(0,*) "KS_GK d:",d," prob:",prob
      
      nK=0 !number of K-stars that are "quiet"
      do 16 i=1,ns
      	if(samp(i).lt.mnoise+merr*sigma)then
      		nK=nK+1
        endif
 16   continue
 
	  write(0,*) "Percentage of Quiet K-stars:",100.0*real(nK)/real(ns)    

      nr=0 !initialize counter to track sample sizes.
      ns=0
      ndAf=0
      ndM=0
      do 17 i=1,nmax
      	if((sptype(i)(1:1).eq.'A').and.(rmag(i).gt.13.0).and.
     .		(rmag(i).lt.14.4))then
      		nr=nr+1
      		ndAf=ndAf+1
      		ref(nr)=noise(i) !store A faint-star noises
      		dataAf(ndAf)=noise(i)
        endif
        if(sptype(i)(1:1).eq.'M')then
        	ns=ns+1
        	ndM=ndM+1
        	samp(ns)=noise(i) !store M star noises
        	dataM(ndM)=noise(i)
        endif
 17   continue
      write(0,*) "Number of A-stars: ",nr, "(faint)"
      write(0,*) "Number of M-stars: ",ns
      
      call kstwo(dataAf,ndAf,dataM,ndM,d,prob)
      write(0,*) "KS_AM d:",d," prob:",prob

C     Calculate mean noise for faint A-stars
      mnoise=0.0
      do 18 i=1,nr
      	mnoise=mnoise+ref(i)
 18   continue
      mnoise=mnoise/real(nr)

C     get the error on the mean from a simple bootstrap    
c      merr=meanboot(seed,ns,samp,nboot,boottemp,nmax,nsamp)
      call avevar(ref,nr,mnoise,merr)
      merr=sqrt(merr)
      write(0,*) "A-star faint noise: ",mnoise,"+/-",merr

      nM=0 !number of K-stars that are "quiet"
      do 19 i=1,ns
      	if(samp(i).lt.mnoise+merr*sigma)then
      		nM=nM+1
        endif
 19   continue
 
	  write(0,*) "Percentage of Quiet M-stars:",100.0*real(nM)/real(ns)
        
      close(10)
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function meanboot(seed,nr,ref,nboot,btemp,nmax,nsamp)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  implicit none
	  integer seed,nr,nboot,nmax,i,j,k
	  real ref(nr),btemp(nboot),nsamp(nmax),ran2,ave,var
	  
	  do 10 i=1,nboot
	  	do 11 k=1,nr
	  		j=ran2(seed)*real(nr)+1.0
	  		nsamp(k)=ref(j) !new sample with replacement
 11		continue
 		btemp(i)=0. !calculate new mean
 		do 12 k=1,nr
 			btemp(i)=btemp(i)+nsamp(k)
 12     continue
 		btemp(i)=btemp(i)/real(nr) !get mean
c 		write(0,*) btemp(i)
 10   continue	  
	 
	  call avevar(btemp,nboot,ave,var)
	  
	  meanboot=sqrt(var)
	 
	  return
	  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n
      REAL ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL s,ep
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      SUBROUTINE kstwo(data1,n1,data2,n2,d,prob)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n1,n2
      REAL d,prob,data1(n1),data2(n2)
C USES probks,sort
      INTEGER j1,j2
      REAL d1,d2,dt,en1,en2,en,fn1,fn2,probks
      call sort(n1,data1)
      call sort(n2,data2)
      en1=n1
      en2=n2
      j1=1 !Next value of data1 to be processed.
      j2=1
      fn1=0.
      fn2=0.
      d=0.
 1    if(j1.le.n1.and.j2.le.n2)then !If we are not done...
         d1=data1(j1)
         d2=data2(j2)
         if(d1.le.d2)then !Next step is in data1.
            fn1=j1/en1
            j1=j1+1
         endif
         if(d2.le.d1)then !Next step is in data2.
            fn2=j2/en2
            j2=j2+1
         endif
         dt=abs(fn2-fn1)
         if(dt.gt.d)d=dt
         goto 1
      endif
      en=sqrt(en1*en2/(en1+en2))
      prob=probks((en+0.12+0.11/en)*d)!Compute significance.
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
      real FUNCTION probks(alam)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      REAL alam,EPS1,EPS2
      PARAMETER (EPS1=0.001, EPS2=1.e-8)
C Kolmogorov-Smirnov probability function.
      INTEGER j
      REAL a2,fac,term,termbf
      a2=-2.*alam**2
      fac=2.
      probks=0.
      termbf=0. !Previous term in sum.
      do 11 j=1,100
         term=fac*exp(a2*j**2)
         probks=probks+term
         if(abs(term).le.EPS1*termbf.or.abs(term).le.EPS2*probks)return
         fac=-fac !Alternating signs in sum.
         termbf=abs(term)
 11   continue
      probks=1. !Get here only by failing to converge.
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE sort(n,arr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
C Sorts an array arr(1:n) into ascending numerical order using the Quicksort algorithm. n
C is input; arr is replaced on output by its sorted rearrangement.
C Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
C auxiliary storage.
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.M)then !Insertion sort when subarray small enough.
         do 12 j=l+1,ir
            a=arr(j)
            do 11 i=j-1,l,-1
               if(arr(i).le.a)goto 2
               arr(i+1)=arr(i)
 11         continue
            i=l-1
 2          arr(i+1)=a
 12      continue
         if(jstack.eq.0)return
         ir=istack(jstack) !Pop stack and begin a new round of partitioning.
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2 !Choose median of left, center, and right elements as partitioning
C                   !element a. Also rearrange so that a(l) ≤ a(l+1) ≤ a(ir).
         temp=arr(k)
         arr(k)=arr(l+1)
         arr(l+1)=temp
         if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l+1).gt.arr(ir))then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l).gt.arr(l+1))then
            temp=arr(l)
            arr(l)=arr(l+1)
            arr(l+1)=temp
         endif
         i=l+1 !Initialize pointers for partitioning.
         j=ir
         a=arr(l+1) !Partitioning element.
 3       continue !Beginning of innermost loop.
            i=i+1 !Scan up to find element > a.
            if(arr(i).lt.a)goto 3
 4          continue
            j=j-1 !Scan down to find element < a.
            if(arr(j).gt.a)goto 4
            if(j.lt.i)goto 5 !Pointers crossed. Exit with partitioning complete.
            temp=arr(i) !Exchange elements.
            arr(i)=arr(j)
            arr(j)=temp
         goto 3 !End of innermost loop.
 5       arr(l+1)=arr(j) !Insert partitioning element.
         arr(j)=a
         jstack=jstack+2
C Push pointers to larger subarray on stack, process smaller subarray immediately.
         if(jstack.gt.NSTACK)pause "NSTACK too small in sort"
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1
      END