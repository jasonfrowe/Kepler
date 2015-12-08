      program sigclipprog
      implicit none
      integer iargc,nmax,npt,nfit,nunit,i,j,npt2
      parameter(nmax=1650000,nfit=18)
      integer tflag(nmax),cflag(nmax)
      double precision sigclip,sol(nfit),serr(nfit,2),Dpvary(nfit),
     .  err(nfit),doe,toff,eoff,phase(nmax),x(nmax),y(nmax),time(nmax),
     .  flux(nmax),ferr(nmax),stdev,std,my,deriv(nmax)
      character*80 filename,fitparsfile,cline
      
      
      if(iargc().lt.2) goto 901
      call getarg(1,filename)
      call getarg(2,cline)
      read(cline,*) sigclip
      if(sigclip.lt.0) goto 901
      if(sigclip.eq.0) sigclip=99.9e30
c      call getarg(2,fitparsfile)
      
      
      nunit=10
      open(unit=nunit,file=filename,status='old',err=902)
      call readkeplc(nunit,nmax,npt,time,flux,ferr)
      close(nunit)
      if(npt.eq.0) goto 999
      
c      nunit=10
c      open(unit=nunit,file=fitparsfile,status='old',err=904)
c      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff,eoff)
c      close(nunit)
c      call marktransit(npt,phase,time,flux,tflag,nfit,sol)

      do 13 i=1,npt
        tflag(i)=0
        cflag(i)=0
 13   continue

      do 12 i=3,13 !upto 10 planets
        if(iargc().ge.i) then
            write(0,*) "Planet ",i-2
            call getarg(i,fitparsfile)
            nunit=10
            open(unit=nunit,file=fitparsfile,status='old',err=904)
           call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff,eoff)
            close(nunit)
            call marktransit(npt,phase,time,flux,tflag,nfit,sol)
        endif
 12   continue

 
      call outlier(npt,time,flux,cflag,deriv)
      
      j=0
      do 10 i=1,npt
        if(tflag(i).eq.0)then
            j=j+1
c            x(j)=time(i)
            y(j)=flux(i)
        endif
 10   continue
      npt2=j
      
      std=stdev(npt2,y,my)
      
      do 11 i=1,npt
        if( (abs(flux(i)-my).lt.sigclip*std).or.(tflag(i).eq.1) )then
            if(cflag(i).eq.0)then
                write(6,*) time(i)-0.5d0+54900.0d0,flux(i),ferr(i)
            endif
        endif
 11   continue
      
      
      goto 999
 901  write(0,*) "Usage: sigclip <filein> <f1.dat> <sigclip>"
      write(0,*) "<filein> Kepler photometry"
      write(0,*) "<sigclip> sigma value for clipping"
      write(0,*) "<f1.dat> Transit solution"
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 904  write(6,*) "Cannot open ",fitparsfile
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine outlier(npt,time,flux,cflag,deriv)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,cflag(npt),i
      double precision time(npt),flux(npt),deriv(npt),std,mean,stdev,
     .  sigcut
      
      sigcut=3.0
      
      do 10 i=1,npt-1
        deriv(i)=flux(i+1)-flux(i)
 10   continue
      
      std=stdev(npt-1,deriv,mean)
      
      do 11 i=1,npt-1
        if((abs(deriv(i)-mean).gt.std*sigcut).and.
     .    (abs(deriv(i+1)-mean).gt.std*sigcut).and.
     .    (deriv(i)/deriv(i+1).lt.0.0d0))then
            cflag(i+1)=1
        endif
 11   continue 
      
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine phasept(npt,time,phase,period,toff)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     npt - number of data points
C     time - times of observations
C     mag - the observations
C     phase - returned phase of observation
C     period - fixed period for data
      implicit none
      integer npt
      double precision time(npt),phase(npt),period,toff

      integer i
      double precision temp

      do 10 i=1,npt
         temp=time(i)
C        Get the phase
         phase(i)=temp/period-int(temp/period)
C        apply optional phase offset to make plot pretty
         phase(i)=phase(i)+toff
C        make sure phase is between 0 and 1
         if(phase(i).lt.0.0) phase(i)=phase(i)+1.0
         if(phase(i).gt.1.0) phase(i)=phase(i)-1.0
 10   continue
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readkeplc(nunit,nmax,npt,dtime,flux,ferr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nmax,npt,i
      double precision dtime(nmax),flux(nmax),ferr(nmax),Keplertime

      Keplertime=54900.0d0

      i=1
      
  9   continue
 10   read(nunit,*,err=9,end=20) dtime(i),flux(i),ferr(i)
         dtime(i)=dtime(i)+0.5d0-Keplertime
c        mag(i)=-2.5*log10(mag(i)+1.0d0)
c        ferr(i)=0.00005
        i=i+1
      goto 10
 20   continue
        
      npt=i-1
      write(0,*)   "-------------------------"
      write(0,500) "Observations read: ",npt
      write(0,*)   "-------------------------"
 500  format(1X,A19,I6)
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nunit,i
      double precision sol(nfit),serr(nfit,2),Dpvary(nfit),toff,
     .  err(nfit),doe
      character*160 command
      
c      write(0,*) "----------------------------"
c      write(0,*) "Reading in starting solution"
c      write(0,*) "----------------------------"
      
      i=0 !initialize line counter
      toff=0.0 !initialize toff to zero
C     Start of loop to read in file
  10  read(nunit,500,end=11) command
 500  format(A160) !command much be contained within first 80 characters
        i=i+1 !increase line counter
        if(command(1:1).eq."#") then
            write(0,*) "Comment line on line ",i
        elseif(command(1:3).eq."SMA") then
            read(command(5:160),*) sol(1),serr(1,1),serr(1,2),Dpvary(1),
     .          err(1)
c           write(0,501) "SMA: ",sol(1),serr(1,1),serr(1,2),Dpvary(1),
c     .          err(1)
        elseif(command(1:3).eq."PMA") then
            read(command(5:160),*) sol(2),serr(2,1),serr(2,2),Dpvary(2),
     .          err(2)
c           write(0,501) "PMA: ",sol(2),serr(2,1),serr(2,2),Dpvary(2),
c     .          err(2)
        elseif(command(1:3).eq."SRA") then
            read(command(5:160),*) sol(3),serr(3,1),serr(3,2),Dpvary(3),
     .          err(3)
c           write(0,501) "SRA: ",sol(3),serr(3,1),serr(3,2),Dpvary(3),
c     .          err(3)
        elseif(command(1:3).eq."PRA") then
            read(command(5:160),*) sol(4),serr(4,1),serr(4,2),Dpvary(4),
     .          err(4)
c           write(0,501) "PRA: ",sol(4),serr(4,1),serr(4,2),Dpvary(4),
c     .          err(4)
        elseif(command(1:3).eq."PER") then
            read(command(5:160),*) sol(5),serr(5,1),serr(5,2),Dpvary(5),
     .          err(5)
c           write(0,501) "PER: ",sol(5),serr(5,1),serr(5,2),Dpvary(5),
c     .          err(5)
        elseif(command(1:3).eq."INC") then
            read(command(5:160),*) sol(6),serr(6,1),serr(6,2),Dpvary(6),
     .          err(6)
c           write(0,501) "INC: ",sol(6),serr(6,1),serr(6,2),Dpvary(6),
c     .          err(6)
        elseif(command(1:3).eq."EPO") then
            read(command(5:160),*) sol(7),serr(7,1),serr(7,2),Dpvary(7),
     .          err(7)
c           write(0,501) "EPO: ",sol(7),serr(7,1),serr(7,2),Dpvary(7),
c     .          err(7)
        elseif(command(1:3).eq."ZPT") then
            read(command(5:160),*) sol(8),serr(8,1),serr(8,2),Dpvary(8),
     .          err(8)
c           write(0,501) "ZPT: ",sol(8),serr(8,1),serr(8,2),Dpvary(8),
c     .          err(8)
        elseif(command(1:3).eq."ALB") then
            read(command(5:160),*) sol(9),serr(9,1),serr(9,2),Dpvary(9),
     .          err(9)
c           write(0,501) "ALB: ",sol(9),serr(9,1),serr(9,2),Dpvary(9),
c     .          err(9)
        elseif(command(1:3).eq."NL1") then
            read(command(5:160),*) sol(10),serr(10,1),serr(10,2),
     .          Dpvary(10),err(10)
c         write(0,501) "NL1: ",sol(10),serr(10,1),serr(10,2),Dpvary(10),
c     .          err(10)
        elseif(command(1:3).eq."NL2") then
            read(command(5:160),*) sol(11),serr(11,1),serr(11,2),
     .          Dpvary(11),err(11)
c         write(0,501) "NL2: ",sol(11),serr(11,1),serr(11,2),Dpvary(11),
c     .          err(11)
        elseif(command(1:3).eq."NL3") then
            read(command(5:160),*) sol(12),serr(12,1),serr(12,2),
     .          Dpvary(12),err(12)
c         write(0,501) "NL3: ",sol(12),serr(12,1),serr(12,2),Dpvary(12),
c     .          err(12)
        elseif(command(1:3).eq."NL4") then
            read(command(5:160),*) sol(13),serr(13,1),serr(13,2),
     .          Dpvary(13),err(13)
c         write(0,501) "NL4: ",sol(13),serr(13,1),serr(13,2),Dpvary(13),
c     .          err(13)
        elseif(command(1:3).eq."ECW") then
            read(command(5:160),*) sol(14),serr(14,1),serr(14,2),
     .          Dpvary(14),err(14)
c         write(0,501) "ECN: ",sol(14),serr(14,1),serr(14,2),Dpvary(14),
c     .          err(14)
        elseif(command(1:3).eq."ESW") then
            read(command(5:160),*) sol(15),serr(15,1),serr(15,2),
     .          Dpvary(15),err(15)
c         write(0,501) "WWW: ",sol(15),serr(15,1),serr(15,2),Dpvary(15),
c     .          err(15)
        elseif(command(1:3).eq."TED") then
            read(command(5:160),*) sol(16),serr(16,1),serr(16,2),
     .          Dpvary(16),err(16)
c         write(0,501) "TED: ",sol(16),serr(16,1),serr(16,2),Dpvary(16),
c     .          err(16)
        elseif(command(1:3).eq."VOF") then
            read(command(5:160),*) sol(17),serr(17,1),serr(17,2),
     .          Dpvary(17),err(17)
        elseif(command(1:3).eq."DIL") then
            read(command(5:160),*) sol(18),serr(18,1),serr(18,2),
     .          Dpvary(18),err(18)
        elseif(command(1:3).eq."OFF") then
            read(command(5:160),*,err=12,end=12) toff,doe
            goto 13
 12         doe=0.0
            read(command(5:160),*) toff
 13     continue
c            write(0,501) "OFF: ",toff,doe
        endif
 501    format(A5,5(1X,PE17.10))
      goto 10 !loop back to read statement
 11   continue !end loop when EOF is reached
c      write(0,*) "----------------------------"
             
      return
      end 
