      program tpsguess
C     Given TPS parameters, guess a transit solution for transitfit
      implicit none
      integer iargc,kid,nunit,i,nmax,npt,kicid,nfit,nper,j,koi
      parameter(nmax=600000,nfit=18)
      real*8 per,mjdepoch,teff,rad,logg,dtime(nmax),mag(nmax),
     .  merr(nmax),itime(nmax),Keplertime,kmag,sol(nfit),M1,M2,Psec,
     .  R1,R2,transitdur,asemi,epoch,boxbin,transitdepth,Dpvary(nfit),
     .  serr(nfit,2),maxtime,toff,err(nfit),doe,incl,tdepth2,td
      character*3 titles(nfit)
      character*80 cline,filename,path,fname
      logical loop
C     Read in physical constants (Pi,G,Msun,..etc.)      
      include "utils/physcons.f"
      
      if(iargc().lt.8) goto 901 !check number of command line arguements
      call getarg(1,cline)
      read(cline,*) kid !read Kelper ID
      call getarg(2,cline)
      read(cline,*) per !read transit period guess
      call getarg(3,cline)
      read(cline,*) mjdepoch !time of first transit (BJD)
      call getarg(4,cline)
      read(cline,*) kmag
      call getarg(5,cline)
      read(cline,*) Teff
      call getarg(6,cline)
      read(cline,*) logg
      call getarg(7,cline)
      read(cline,*) rad
      call getarg(8,cline)
      read(cline,*) td
      koi=-1
      if(iargc().ge.10)then
        call getarg(10,cline)
        read(cline,*) koi
      endif
      
      nper=1 !okay to fit period
      if(iargc().ge.9)then
        call getarg(9,cline)
        read(cline,*) nper
      endif
      
     
c      Keplertime=54900.0
cc  add HJD 5/2/09 adjustment for field center
c      epoch=mjdepoch+0.50020-Keplertime
cc  add time dependent change for field center
c      epoch=epoch+4.1d-5*(epoch-53.0)

C     BJD Feb 5,2010 correction
c      Keplertime=54900.0d0
      Keplertime=67.0
      epoch=mjdepoch-Keplertime
c      write(0,*) epoch,mjdepoch,Keplertime
c      read(5,*)
c      epoch=epoch+sin(2.0d0*pi*epoch/372.5d0-1.1208d0)*0.00280758d0
c      epoch=epoch+0.5d0 !MJD half day offset
c      epoch=mjdepoch

C     Get filename for Kepler Q1 lightcurve      
      if(kID.lt.10)then
        write(filename,507) "klc0000000",kID,".d.dat"
 507    format(A10,I1,A6)      
      elseif(kID.lt.100)then
        write(filename,506) "klc000000",kID,".d.dat"
 506    format(A9,I2,A6)      
      elseif(kID.lt.1000)then
        write(filename,505) "klc00000",kID,".d.dat"
 505    format(A8,I3,A6)      
      elseif(kID.lt.10000)then
        write(filename,504) "klc0000",kID,".d.dat"
 504    format(A7,I4,A6)      
      elseif(kID.lt.100000)then
        write(filename,503) "klc000",kID,".d.dat"
 503    format(A6,I5,A6)      
      elseif(kID.lt.1000000)then
        write(filename,502) "klc00",kID,".d.dat"
 502    format(A5,I6,A6)
      elseif(kID.lt.10000000)then
        write(filename,501) "klc0",kID,".d.dat"
 501    format(A4,I7,A6)
      else
        write(filename,500) "klc",kID,".d.dat"
 500    format(A3,I8,A6)
      endif
      
c      path="/media/Streetsville/rowe/Kepler/q3phot/" !where the files are
c      path="/media/Etobicoke/Kepler/q3phot/"
      path="./"
      i=2
c      i=80
c      loop=.true.
c      do while(loop)
c        if(path(i:i).eq." ") then
c            i=i-1
c        else
c            loop=.false.
c        endif
c      enddo
      
      nunit=10  !read in photometry
c      open(unit=nunit,file=path(1:i)//filename(1:15),status='old',
c     .  err=902)
      open(unit=nunit,file=filename,status='old',
     .  err=902)
      call readkeplc(nunit,nmax,npt,dtime,mag,merr,itime,Keplertime)
      close(nunit)
      
      maxtime=dtime(1)
      do 14 i=2,npt
        maxtime=max(maxtime,dtime(i))
 14   continue

C     apply a boxcar filter
c      boxbin=2.0 !filter width (days)
c      call boxcar(npt,dtime,mag,merr,boxbin)
      
c      filename="/media/Streetsville/rowe/Kepler/tpsq3/kic_0_16.5.txt" !read in KIC information
c      filename="/media/Etobicoke/Kepler/tpsq3/kic_0_16.5.txt"
c      filename="/home/rowe/p555/transitfit/Kepler/kic_0_16.5.txt"
      filename="kic_0_16.5.txt"
c      open(unit=nunit,file=filename,status='old',err=902)
c
c      i=0
c      j=1
c 10   read(nunit,*,end=11,err=903) kicid,kmag,teff,logg,rad
c        j=j+1
c        if(kid.eq.kicid)then
c            i=1
c            goto 11
c        endif
c        goto 10
c 11   continue
c      if(i.eq.0)then
c        kmag=0.0
c        teff=0.0
c        rad=0.0
c        logg=0.0
c      endif
c      if(kmag.eq.0) kmag=0.0
c      if(rad.eq.0) rad=1.0
c      if(logg.eq.0.0) logg=4.5
c      close(nunit)
      
c      write(0,508) kid,kmag,int(teff),logg,rad
 508  format(I8,1X,F6.3,1X,I5,1X,F6.3,1X,F5.3)
 
C     Start building solution
            
C     Get mass from logg
      sol(1)=10.0**logg*(rad*rad*Rsun*Rsun)/(Msun*G)/100.0
c      write(0,*) "sol(1):",sol(1)
C     Assume a mass of 0 Jupiter
      sol(2)=0.0d0
c      write(0,*) "sol(2):",sol(2)
C     Radius from KIC
      sol(3)=rad
c      write(0,*) "sol(3):",sol(3)
C     Period from TPS
      sol(5)=per
      
      M1=sol(1)*Msun !kg ; mass of star
      M2=sol(2)*Mjup !kg ; mass of planet
      Psec=sol(5)*24.0*60.0*60.0 !sec ; period of planet
      asemi=(Psec*Psec*G*(M1+M2)/(4.0*Pi*Pi))**(1.0/3.0) !m
      R1=sol(3)*Rsun  !radius of star
      R2=0.0!assume zero for now. sol(4)*Rjup  !radius of planet     
      transitdur=Psec*(R1+R2)/(Pi*asemi)/(60.0*60.0*24.0) !in days
c      write(6,*) "Transit Duration (h):",transitdur*24.0
      incl=acos(0.5d0*R1/asemi)*180.0d0/Pi

      write(0,*) "td:",td
      if(td.le.1.0e-6) td=tdepth2(npt,dtime,mag,transitdur,epoch,per)
      write(0,*) "td:",td
      if(td.gt.0.0d0)then
        sol(4)=sol(3)*Rsun/Rjup*sqrt(td)!*2.0/3.0
      else
        sol(4)=0.02
      endif
        sol(4)=max(abs(sol(4)),0.02)
c      write(0,*) "sol(4):",sol(4) !write out planet radius
c      write(0,*) "sol(5):",sol(5) !write out planet period
      sol(6)=incl !incliation
c      write(0,*) "sol(6):",sol(6)
      toff=0.75-(epoch/per-int(epoch/per))
      if(toff.lt.0.0)toff=toff+1.0
C     **Change from Phi to Epoch
c      sol(7)=pi-2.0*pi*(epoch/per-int(epoch/per))
c      if(sol(7).lt.0.0)sol(7)=sol(7)+2.0*pi
      sol(7)=epoch
c      write(0,*) "sol(7):",sol(7)
      sol(8)=0.0d0
      sol(9)=0.0d0
      sol(10)= 0.410769d0
      sol(11)=-0.108909d0
      sol(12)= 0.904020d0
      sol(13)=-0.437364d0
      sol(14)=0.0
      sol(15)=0.0
      sol(16)=0.0
      sol(17)=0.0
      sol(18)=0.0
      
      Dpvary(1)=0.02*sol(1)
      Dpvary(2)=0.02*sol(2)
      Dpvary(3)=0.1*sol(3)
      Dpvary(4)=0.1*sol(4)
      Dpvary(5)=4.0d-3
      Dpvary(6)=0.5
      Dpvary(7)=0.02
      Dpvary(8)=1.0d-5
      Dpvary(9)=0.2
      Dpvary(10)=0.005
      Dpvary(11)=0.005
      Dpvary(12)=0.005
      Dpvary(13)=0.005
      Dpvary(14)=0.05
      Dpvary(15)=10.0
      Dpvary(16)=10.0
      Dpvary(17)=1.0
      Dpvary(18)=0.01
      
C     Default is to disable priors
      do 13 i=1,nfit
        serr(i,1)= 0.0d0
        serr(i,2)= 0.0d0
        err(i)=0.0
 13   continue
      serr(4,2)=-1.0
      if(epoch+per.lt.maxtime) serr(5,2)=-1.0
      if(nper.eq.0) serr(5,2)=0.0d0
      serr(6,2)=-1.0
      serr(7,2)=-1.0
      serr(8,2)=-1.0
c      serr(16,2)=-1.0

      doe=0.0
      
      call exportfit(nfit,sol,serr,Dpvary,err,doe,toff,titles)
      
      open(unit=11,file="plotend.dat")
      write(11,*) kmag
      write(11,*) int(Teff)
      write(11,*) logg
      write(11,*) rad
      close(11)
      
C     Get filename for Kepler Q1 lightcurve      
      if(kID.lt.10)then
        write(filename,510) "klc0000000",kID,".d.dat"
 510    format(A10,I1,A6)      
      elseif(kID.lt.100)then
        write(filename,511) "klc000000",kID,".d.dat"
 511    format(A9,I2,A6)      
      elseif(kID.lt.1000)then
        write(filename,512) "klc00000",kID,".d.dat"
 512    format(A8,I3,A6)      
      elseif(kID.lt.10000)then
        write(filename,513) "klc0000",kID,".d.dat"
 513    format(A7,I4,A6)      
      elseif(kID.lt.100000)then
        write(filename,514) "klc000",kID,".d.dat"
 514    format(A6,I5,A6)      
      elseif(kID.lt.1000000)then
        write(filename,515) "klc00",kID,".d.dat"
 515    format(A5,I6,A6)
      elseif(kID.lt.10000000)then
        write(filename,516) "klc0",kID,".d.dat"
 516    format(A4,I7,A6)
      else
        write(filename,517) "klc",kID,".d.dat"
 517    format(A3,I8,A6)
      endif
      
C     Get filename for Kepler Q1 lightcurve      
      if(kID.lt.10)then
        write(fname,519) "f_",kID,".dat"
 519    format(A2,I1,A4)      
      elseif(kID.lt.100)then
        write(fname,520) "f_",kID,".dat"
 520    format(A2,I2,A4)      
      elseif(kID.lt.1000)then
        write(fname,521) "f_",kID,".dat"
 521    format(A2,I3,A4)      
      elseif(kID.lt.10000)then
        write(fname,522) "f_",kID,".dat"
 522    format(A2,I4,A4)      
      elseif(kID.lt.100000)then
        write(fname,523) "f_",kID,".dat"
 523    format(A2,I5,A4)      
      elseif(kID.lt.1000000)then
        write(fname,524) "f_",kID,".dat"
 524    format(A2,I6,A4)
      elseif(kID.lt.10000000)then
        write(fname,525) "f_",kID,".dat"
 525    format(A2,I7,A4)
      else
        write(fname,526) "f_",kID,".dat"
 526    format(A2,I8,A4)
      endif      
      
      open(unit=11,file="plot.dat")
      write(11,518) filename
 518  format(A17)
      write(11,*) "1"
      write(11,527) fname
 527  format(A14)
      write(11,*) koi
      write(11,*) kid
      write(11,*) kmag
      write(11,*) int(Teff)
      write(11,*) logg
      write(11,*) rad
      write(11,*) "0.0"
      write(11,*) "1.0"
      write(11,*) "0.0"
      close(11)
      
      goto 999
 901  write(0,*) "Usage: tpsguess <kID> <Period> <Epoch>,kmag,T,lg,R,Td"
      write(0,*) "  "
      write(0,*) "<kID>: Kepler ID"
      write(0,*) "<Period>: Periodicity of Transits (days)"
      write(0,*) "<Epoch>: Time of first transit (MJD)"
      goto 999
 902  write(0,*) "Cannot open",filename
      goto 999
 903  write(0,*) "Error in TPS line: ",j
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function tdepth2(npt,time,mag,tdur,epoch,per)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,j
      double precision time(npt),mag(npt),tdur,epoch,per,toff,ph,mean,
     .  ph1,ph2
     
      toff=0.75-(epoch/per-int(epoch/per))
      ph1=0.75-0.5d0*tdur/per
      if(ph1.lt.0.5)ph1=0.5
      ph2=0.75+0.5d0*tdur/per
      if(ph2.gt.1.0)ph2=1.0
c      write(0,*) epoch,per
      
      mean=0.0d0
      j=0
      do 10 i=1,npt
        ph=(time(i)/per-int(time(i)/per))+toff
        if(ph.lt.0.0) ph=ph+1.0
        if(ph.gt.1.0) ph=ph-1.0
c        write(0,*) time(i),ph,mag(i)
        if((ph.ge.ph1).and.(ph.le.ph2))then
            j=j+1
            mean=mean+mag(i)
c            write(0,*) ph,mag(i)
        endif
 10   continue
      if(j.gt.0) mean=mean/dble(j)
c      if(mean.le.0) mean=1.0e-5
      write(0,*) "mean:",mean,ph1,ph2
      tdepth2=abs(1.0-10.0**(mean/-2.5d0))
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function transitdepth(npt,time,mag,transitdur,epoch)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i1,i2,i
      real*8 time(npt),mag(npt),transitdur,epoch,td2,fmax
      logical loop
      
      td2=transitdur/2.0d0 !half of the transit duration
      i1=1
      loop=.true.
      do while(loop)
        if(epoch-time(i1).lt.td2)then
            i1=i1-1
            loop=.false.
        else
            i1=i1+1
        endif
        if(i1.gt.npt)then
            loop=.false.
            i1=npt
        endif
      enddo
      if(i1.lt.1) i1=1
      
      i2=1
      loop=.true.
      do while(loop)
        if(time(i2)-epoch.gt.td2)then
            loop=.false.
        else
            i2=i2+1
        endif
        if(i2.gt.npt)then
            loop=.false.
            i2=npt
        endif
      enddo
      
c      write(0,*) i1,time(i1),epoch
c      write(0,*) i2,time(i2),epoch
      
      fmax=0.0
      do 10 i=i1,i2
        fmax=max(fmax,mag(i))
 10   continue
      transitdepth=1.0-10.0**(fmax/-2.5d0)
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mag2flux(npt,mag,flux)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      double precision mag(npt),flux(npt)
      
      do 10 i=1,npt
        flux(i)=10.0**(mag(i)/-2.5)
 10   continue
 
      return
      end
