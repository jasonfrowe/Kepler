      program corotplot
      implicit none
      integer nmax,np,i,j,k
      parameter(nmax=21210,np=5)
      integer npt(np),id(np)
      real rmag(np),teff(np),rms(3,np),expt(3,np),tmin,tmax,mmin,mmax,
     .	xp(nmax),yp(nmax)
      double precision time(nmax,np),mag(nmax,np),merr(nmax,np),
     .	tt(nmax),tm(nmax),te(nmax)
      character*3 sp(np),cl(np)
      character*80 filename,corotfile
      
      tmin=2691.75312436
      tmax=2843.76648217
      
      filename="corot_12rms.dat"
      
      open(unit=10,file=filename,status='old',err=901)
      
      call pgopen('?')
      call pgsch(3.0)
      call pgsubp(1,np)
      call pgvport(0.05,0.85,0.1,0.85)
      
      i=1
      k=0
 10   read(10,*,end=11) id(i),rmag(i),teff(i),sp(i),cl(i),
     .	(rms(j,i),j=1,3),(expt(j,i),j=1,3)
      	write(corotfile,501) id(i),".asc"
 501  	format(I5,A4)
 
      	call readcorotdata(corotfile,nmax,np,i,npt,time,mag,merr,tt,tm,
     .		te)
      	if(i.eq.np)then
      		call findminmax(nmax,np,i,npt,mag,mmin,mmax)
 			call plotcorot(nmax,np,i,npt,time,mag,merr,tmin,tmax,mmin,
     .			mmax,xp,yp,id,rmag,teff,sp,cl,rms,expt)
      	endif
        call pgpage()
      
      	i=i+1
      	if(i.gt.np) i=1
      goto 10
 11   continue
      close(10)
      i=i-1
      call findminmax(nmax,np,i,npt,mag,mmin,mmax)
 	  call plotcorot(nmax,np,i,npt,time,mag,merr,tmin,tmax,mmin,
     .	mmax,xp,yp,id,rmag,teff,sp,cl,rms,expt)
      
      call pgclos()
      
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  subroutine plotcorot(nmax,np,npr,npt,time,mag,merr,tmin,tmax,mmin,
     .	mmax,xp,yp,id,rmag,teff,sp,cl,rms,expt)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nmax,np,npr,npt(np),i,j,id(np)
      real xp(nmax),yp(nmax),tmin,tmax,mmin,mmax,rmag(np),teff(np),
     .	rms(3,np),expt(3,np)
      double precision time(nmax,np),mag(nmax,np),merr(nmax,np)
      character*80 line1,line2,line3
      character*3 sp(np),cl(np)

      do 10 i=1,npr
      	call pgpanl(1,i)
      	call pgwindow(tmin,tmax,mmax,mmin)
      	call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
        if(i.eq.np/2+1) call pglabel("","mmag","")
        if(i.eq.npr) call pglabel("Days","","")
      	do 11 j=1,npt(i)
      		xp(j)=real(time(j,i))
      		yp(j)=real(mag(j,i))
 11		continue
 	  	call pgpt(npt(i),xp,yp,-1)
c 	  	write(line1,500) id(i),rmag(i)
c 500    format(I5,1X,F6.3)
        write(line1,503) "Rmag= ",rmag(i)
 503    format(A6,F6.3)
 	    call pgptxt(tmax+0.01*(tmax-tmin),mmin+0.16*(mmax-mmin),
     .		0.0,0.0,line1)
        write(line2,501) teff(i),sp(i),cl(i)
 501    format(F5.0,1X,A3,1X,A3)
        call pgptxt(tmax+0.01*(tmax-tmin),mmin+0.32*(mmax-mmin),
     .		0.0,0.0,line2)
        write(line3,502) rms(1,i),expt(1,i)
 502    format(F10.4,1X,F10.4)
        call pgptxt(tmax+0.01*(tmax-tmin),mmin+0.48*(mmax-mmin),
     .		0.0,0.0,line3)
        write(line3,502) rms(2,i),expt(2,i)
        call pgptxt(tmax+0.01*(tmax-tmin),mmin+0.64*(mmax-mmin),
     .		0.0,0.0,line3)
        write(line3,502) rms(3,i),expt(3,i)
        call pgptxt(tmax+0.01*(tmax-tmin),mmin+0.8*(mmax-mmin),
     .		0.0,0.0,line3)
 10   continue
 
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findminmax(nmax,np,npr,npt,mag,mmin,mmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nmax,np,npt(np),i,npr,j
      double precision mag(nmax,np)
      real mmin,mmax
      
      mmin=real(mag(1,1))
      mmax=real(mag(1,1))
      do 11 i=1,npr
      	do 10 j=2,npt(i)
      		mmin=min(mag(j,i),mmin)
      		mmax=max(mag(j,i),mmax)
 10   	continue
 11   continue
 
      return
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readcorotdata(corotfile,nmax,np,i,npt,time,mag,merr,
     .	tt,tm,te)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nmax,np,i,npt(np),j,nbin
      double precision time(nmax,np),mag(nmax,np),merr(nmax,np),tbin,
     .	tt(nmax),tm(nmax),te(nmax)
      character*80 corotfile
      character dumc
      
      open(unit=11,file=corotfile,status='old',err=901)
      read(11,*) dumc !first line of line is a comment
      
      j=1
  10  read(11,*,end=11) time(j,i),mag(j,i),merr(j,i)
        mag(j,i)=mag(j,i)*1000.0
        j=j+1
      goto 10
  11  continue
      npt(i)=j-1
      
c      tbin=6.0*60.0
c      nbin=npt(i)
c      do 12 j=1,nbin
c      	tt(j)=time(j,i)
c      	tm(j)=mag(j,i)
c      	te(j)=merr(j,i)
c 12   continue
c      call bindt(nbin,tt,tm,te,tbin)
c      npt(i)=nbin
c      do 13 j=1,npt(i)
c      	time(j,i)=tt(j)
c      	mag(j,i)=tm(j)
c      	merr(j,i)=te(j)
c 13   continue
      
      goto 999
 901  write(6,*) "Cannot open ",corotfile
      pause
      goto 999
 999  return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindt(npt,time,mag,merr,tbin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,nbins,bin,nmax,j
      parameter(nmax=650000)
      real*8 time(npt),mag(npt),merr(npt),tbin,tmin,tmax,ltime,
     .     avgm(nmax),avgt(nmax),avge(nmax),abin(nmax),stdev(nmax),
     .     var(nmax),ep(nmax),s,p,avgs(nmax),sigcut


      sigcut=1.0

C     assume time is in days and tbin is in minutes

      do 5 i=1,nmax
         avgm(i)=0.0
         avge(i)=0.0
         avgt(i)=0.0
         abin(i)=0.0
 5    continue

      tmin= 99.9e30
      tmax=-99.9e30
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