CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindt(npt,time,mag,merr,tbin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,nbins,bin,nmax,j
      parameter(nmax=650000)
      double precision time(npt),mag(npt),merr(npt),tbin,tmin,tmax,
     .  ltime,avgm(nmax),avgt(nmax),avge(nmax),abin(nmax),stdev(nmax),
     .  var(nmax),ep(nmax),s,p,avgs(nmax),sigcut


      sigcut=1.0

C     assume time is in days and tbin is in minutes

      do 5 i=1,nmax
         avgm(i)=0.0
         avge(i)=0.0
         avgt(i)=0.0
         abin(i)=0.0
 5    continue

      tmin=time(1)
      tmax=time(1)
      do 10 i=1,npt
         tmin=min(tmin,time(i))
         tmax=max(tmax,time(i))
 10   continue

      ltime=tmax-tmin
      nbins=int(ltime*24.0*60.0/tbin+0.5)+1

      do 20 i=1,npt
         bin=int(dble(nbins)*(time(i)-tmin)/ltime)+1
         if(bin.gt.nmax) write(6,*) "WARNING, nmax too small in bindt"
         avgm(bin)=avgm(bin)+mag(i)/merr(i)
         avgt(bin)=avgt(bin)+time(i)/merr(i)
         avge(bin)=avge(bin)+1.0
         abin(bin)=abin(bin)+1.0/merr(i)
 20   continue

      do 21 i=1,nbins
         avgs(i)=avgm(i)/abin(i)
         ep(i)=0.
         var(i)=0.
         abin(i)=0
 21   continue

      do 22 i=1,npt
         bin=int(dble(nbins)*(time(i)-tmin)/ltime)+1
         s=mag(i)-avgs(bin)
         ep(bin)=ep(bin)+s
         p=s*s
         var(bin)=var(bin)+p
         abin(bin)=abin(bin)+1.0
 22   continue

      do 23 i=1,nbins
         var(i)=(var(i)-ep(i)**2/abin(i))/(abin(i)-1)
         stdev(i)=sqrt(var(i))
         if(abin(i).eq.0) stdev(i)=0.
c         write(6,*) avgs(i),stdev(i)
 23   continue

      do 6 i=1,nmax
         avgm(i)=0.0
         avge(i)=0.0
         avgt(i)=0.0
         abin(i)=0.0
 6    continue

      do 24 i=1,npt
         bin=int(dble(nbins)*(time(i)-tmin)/ltime)+1
         if(abs(mag(i)-avgs(bin)).lt.sigcut*stdev(bin)) then
            avgm(bin)=avgm(bin)+mag(i)/merr(i)
            avgt(bin)=avgt(bin)+time(i)/merr(i)
            avge(bin)=avge(bin)+1.0
            abin(bin)=abin(bin)+1.0/merr(i)
         endif
 24   continue

      j=0
      do 30 i=1,nbins
         if(abin(i).gt.0.0) then
            j=j+1
            avgm(i)=avgm(i)/abin(i)
            avgt(i)=avgt(i)/abin(i)
            avge(i)=(avge(i)**0.5)/abin(i)
            time(j)=avgt(i)
            mag(j)=avgm(i)
            merr(j)=avge(i)
c            write(6,*) time(j),mag(j),merr(j)
         endif
 30   continue
      npt=j
      
      return
      end