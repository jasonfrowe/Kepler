CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readfreqs(nunit,ma,nstar,a,sns,ztime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,i,ma,nstar,k,j,cc,nmax
      parameter(nmax=650000)
      double precision ztime,a(nmax),twopi,pi,sns(nmax)

      pi=3.141592654
      twopi=2.0*pi

      read(nunit,*) ztime
      read(nunit,*) ma,nstar
      j=(ma-nstar)/nstar
      do 10 k=1,nstar
         cc=(j+1)*(k-1)+2
         read(nunit,*) a(cc-1),(a(i),i=cc,cc+j-2,2),
     .        (a(i+1),i=cc,cc+j-2,2),sns(k)
c         write(6,500) a(cc-1),(a(i),i=cc,cc+j-2,2),
c     .        (a(i+1),i=cc,cc+j-2,2),sns(k)
         a(cc-1)=a(cc-1)*twopi
 500     format(20(E14.7,1X))
 10   continue

      return
      end
