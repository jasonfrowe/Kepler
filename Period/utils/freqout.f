CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine freqout(ma,nstar,a,aerr,sns)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ma,nstar,nfmax,i,j,k,cc,b,nstarmax
      parameter(nfmax=2000,nstarmax=300)
      double precision a(nfmax),twopi,pi,aerr(nfmax),sns(nstarmax)

      pi=3.141592654
      twopi=2.0*pi

      j=(ma-nstar)/nstar
      do 38 k=1,nstar
         cc=(j+1)*(k-1)+2
         do 39 i=cc,cc+j-2,2
            if(a(i).lt.0.0) then
               a(i)=abs(a(i))
               a(i+1)=a(i+1)+pi
            endif
            b = int(a(i+1)/twopi)
            a(i+1)=a(i+1)-b*twopi
            if (a(i+1).lt.0) a(i+1)=twopi+a(i+1)
 39      continue
 38   continue

      write(14,502) ma,nstar
      
      j=(ma-nstar)/nstar
      do 10 k=1,nstar
         cc=(j+1)*(k-1)+2
         write(6,504) a(cc-1)/twopi,(a(i),i=cc,cc+j-2,2),
     .        (a(i+1),i=cc,cc+j-2,2),sns(k)
         write(14,504) a(cc-1)/twopi,(a(i),i=cc,cc+j-2,2),
     .        (a(i+1),i=cc,cc+j-2,2),sns(k)
         write(6,504) aerr(cc-1),(aerr(i),i=cc,cc+j-2,2),
     .        (aerr(i+1),i=cc,cc+j-2,2)
 10   continue

 500  format(20(E14.7,1X))
c 500  format(20(F13.6,1X))
 501  format(A3,1X,F11.4)
 502  format(I4,1X,I4)
 504  format(17(1PE17.10,1X))

      return
      end
