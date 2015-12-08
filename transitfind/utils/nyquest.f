CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function nyquest(nmax,npt,time,dt,p)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,ndt,i,p(nmax)
      double precision time(npt),dt(nmax),mindt,maxdt,mode,modecal,mean,
     .  std,stdev,sigcut

      sigcut=1.0

      ndt=0
      mindt= 99.9e30
      maxdt=-99.9e30
      mean=0.
      do 10 i=2,npt
         ndt=ndt+1
         dt(ndt)=time(i)-time(i-1)
         mindt=min(dt(ndt),mindt)
         maxdt=max(dt(ndt),maxdt)
         mean=mean+dt(ndt)
 10   continue
      mean=mean/dble(ndt)
c      std=stdev(ndt,dt,mean)
c      mindt=max(mindt,mean-sigcut*std)
c      maxdt=min(maxdt,mean+sigcut*std)

      call rqsort(npt,dt,p)
      mode=dt(p(npt/2))

c      mode=modecal(ndt,dt,mindt,maxdt)
      
c      write(6,*) "mean,mode,mindt,maxdt:",mean,mode,mindt,maxdt

      nyquest=1.0/(2.0*mode)

      return
      end