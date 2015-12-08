CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fixa(ma,nstar,a)
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


      return
      end
