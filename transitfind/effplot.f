      program effplot
      implicit none
      integer nmax,iargc,nunit,nsize,nPer,i,j,k,ntotal,nfind
      parameter(nmax=500000,nsize=5,nPer=600)
      integer npt(nsize),samples(nPer,nsize),passes(nPer,nsize)
      real psize(nsize),x(nPer),y(nPer),Pin,Pout,PindPout,SN,Rp,Epo,ph,
     .  cut,Epoout,epcheck,Pmax,z(nPer),Pmin
      character*80 filename
      data psize/0.06,0.09,0.18,0.45,0.89/
      
      ntotal=0
      nfind=0
      
      cut=0.01
      
      do 15 i=1,nPer
        z(i)=1.0
        do 16 j=1,nsize
            samples(i,j)=0
            passes(i,j)=0
 16     continue
 15   continue
      
      if(iargc().lt.1) goto 901 !check number of command line arguements
      call getarg(1,filename) !get filename for input data
      nunit=10
      
      open(unit=nunit,file=filename,status='old',err=902)
      
      Pmax=2.0
      Pmin=7000.0
 10   read(nunit,*,end=11) Pin,Pout,PindPout,SN,Rp,Epo,Epoout  
        Pmax=max(Pmax,Pin) !find maximum period scanned 
        Pmin=min(Pmin,Pin)
        if(PindPout.lt.0.9) PindPout=1.0/PindPout
        i=0
        do 12 i=1,nsize
            if(Rp.eq.psize(i)) j=i
 12     continue
        if(i.gt.0)then
            samples(int(Pin),j)=samples(int(Pin),j)+1
            ntotal=ntotal+1 !count total trails
            ph=PindPout-int(PindPout)
            epcheck=(Epo-Epoout)/Pout-int((Epo-Epoout)/Pout)
            if(epcheck.lt.0.0) epcheck=epcheck+1.0
            if(epcheck.gt.0.25) epcheck=epcheck-0.5 !catch half periods
            if(epcheck.gt.0.25) epcheck=epcheck-0.5 !catch whole periods
            if(ph.gt.0.5) ph=ph-1.0
            if((abs(ph).lt.cut).and.(PindPout.lt.25.0).and.
     .        (abs(epcheck).lt.0.02))then
                passes(int(Pin),j)=passes(int(Pin),j)+1
                nfind=nfind+1 !count total trails that pass
            endif
        endif
      goto 10
 11   continue

      close(nunit)
      
      
      call pgopen('?')
      call PGPAP ( 6.0 ,1.0)  
c      call PGENV ( -0.5 , 0.5 , -0.5 , 0.5 , 1 , -2 )  
      call pgwindow(Pmin,Pmax,0.0,1.0)
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Period (days)","Completeness","")
      
      do 13 i=1,nsize
        call pgsci(i)
        k=0
        do 14 j=1,nPer
            if(samples(j,i).gt.0)then
                k=k+1
                x(k)=real(j)
                y(k)=real(passes(j,i))/real(samples(j,i))
c                write(0,*) k,x(k),y(k)
            endif
 14     continue
c        call pgline(k,x,y)
        call pgpt(k,x,y,4)
        call bindt(k,x,y,z,10.0)
        call pgline(k,x,y)       
 13   continue
 
      call pgclos()

      write(0,*) "global counts",nfind,ntotal
      write(0,*) "global rate: ",real(nfind)/real(ntotal)
      
      goto 999
 901  write(0,*) "Usage: effplot eff.dat"
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindt(npt,time,mag,merr,tbin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,nbins,bin,nmax,j
      parameter(nmax=650000)
      real time(npt),mag(npt),merr(npt),tbin,tmin,tmax,
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
c      nbins=int(ltime*24.0*60.0/tbin+0.5)+1
      nbins=int(ltime/tbin+0.5)+1

      do 20 i=1,npt
         bin=int(real(nbins)*(time(i)-tmin)/ltime)+1
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
         bin=int(real(nbins)*(time(i)-tmin)/ltime)+1
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
         bin=int(real(nbins)*(time(i)-tmin)/ltime)+1
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