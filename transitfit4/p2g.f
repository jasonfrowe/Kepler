      program phys2geo
      implicit none
C     Convert a Physical space model to a geometric model
      integer nfitg,nfit,nunit,iargc,i
      parameter(nfitg=18,nfit=18)
      double precision solg(nfitg),serrg(nfitg,2),Dpvaryg(nfitg),
     .  errg(nfitg),doe,toff,sol(nfit),serr(nfit,2),err(nfit),Psec,
     .  M1,M2,R1,aConst,asemi,Pi,Pid2,G,sb,Msun,Mearth,Mjup,Rjup,Rsun,
     .  Lsun,incl,R2,eccn,tPi,K,AU
      character*3 titles(nfit)
      character*80 inputsol      

      if(iargc().lt.1) goto 901
      
C     Constants
      Pi=acos(-1.d0)   !Pi
      Pid2=Pi/2.0d0
      tPi=Pi*2.0d0
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      sb=5.6704d-8 !W m^-2 K^-4 Stefan-Boltzman constant
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Mjup=317.833d0*Mearth !kg  mass of Jupiter  
      Rjup=142980.0d0*1000.0d0/2.0d0 !m  radius of Jupiter     
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Lsun=3.839d26 !W Solar Luminosity  
      AU=1.49598e11    
      
      call getarg(1,inputsol)
      nunit=10
      open(unit=nunit,file=inputsol,status='old',err=901)
      call getgfitpars(nunit,nfitg,solg,serrg,Dpvaryg,errg,doe,toff)
      
      do 10 i=1,nfit
        serr(i,1)=0.0d0
        serr(i,2)=0.0d0
        err(i)=0.0d0
 10   continue
      do 11 i=1,6
        serr(i,2)=-1.0d0
 11   continue
      serr(15,2)=-1.0d0
      
      Psec=solg(5)*8.64d4 !sec ; period of planet
      M1=solg(1)*Msun !kg ; mass of star
      M2=solg(2)*Mjup !kg ; mass of planet
      R1=solg(3)*Rsun !m ; stellar radius
      R2=solg(4)*Rjup !m ; planet radius
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      asemi=(M1+M2)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst
c      write(6,*) "asemi: ",asemi/AU,asemi/R1
      incl=solg(6)
      if(incl.gt.90.0d0)incl=180.0-incl     
      incl=Pi*(incl)/180.0d0
      
      eccn=sqrt(solg(14)*solg(14)+solg(15)*solg(15)) !eccentricity
      if(eccn.ge.1.0) eccn=0.99
c      if(eccn.eq.0.0d0)then
c        w=0.0d0
c      else
c        w=atan(sol(8)/sol(7))
c        if((sol(14).gt.0.0d0).and.(sol(15).lt.0.0d0))then
c            w=tPi+w
c        elseif((sol(14).lt.0.0d0).and.(sol(15).gt.0.0d0))then 
c            w=Pi+w
c        elseif((sol(14).lt.0.0d0).and.(sol(15).lt.0.0d0))then
c            w=Pi+w
c        endif
c      endif

      K=2.0*pi*G*M2**3*(sin(incl))**3/
     .  (Psec*(1.0d0-eccn*eccn)**(3.0d0/2.0d0)*(M1+M2)*(M1+M2))
      K=K**(1.0d0/3.0d0)
c      write(6,*) "K: ",K
      
      sol(1)=solg(7)  !Epoch
c      write(6,*) sol(1)
      sol(2)=solg(5)  !Period
c      write(6,*) sol(2)
      sol(3)=(asemi*cos(incl)/R1)**2.0 !b^2 
      if(sol(3).lt.0.1)then
c        write(6,*) 0.1
      else
c        write(6,*) sol(3)
      endif 
      sol(4)=R2/R1  !R/R*
c      write(6,*) sol(4)
      sol(5)=asemi/R1*tPi/solg(5)/sqrt( (1+sol(4))**2.0d0-sol(3))/
     .  (1-solg(15))*sqrt(1-eccn*eccn) !z/R*
c      sol(5)=asemi/R1*tPi/solg(5)/sqrt(1-sol(3))/(1-solg(15))
c      write(6,*) sol(5)
      sol(6)=solg(8) !zpt
c      write(6,*) sol(6)
      sol(7)=solg(14) !ecw
c      write(6,*) sol(7)
      sol(8)=solg(15) !esw
c      write(6,*) sol(8)
      sol(9)=K !radial velocity amplitude
c      write(6,*) sol(9)
      sol(10)=solg(17) !Velocity zero point
c      write(6,*) sol(10)
      sol(11)=solg(10)
      sol(12)=solg(11)
      sol(13)=solg(12)
      sol(14)=solg(13)
c      write(6,*) sol(11),sol(12),sol(13),sol(14)
      sol(15)=solg(16) !TED
c      write(6,*) sol(15)
      sol(16)=0.0 !ellipsoidal variations 
      sol(17)=0.0 !Phase changes of planet (Albedo)
      sol(18)=solg(18) !Dilution
c      write(6,*) sol(16)
c      write(6,*) sol(17)
c      write(6,*) sol(18)
      
      if(sol(3).lt.0.1) sol(3)=0.1
      
      call exportfit(nfit,sol,serr,err,titles)
      
      goto 999
 901  write(0,*) "Cannot open" ,inputsol
      goto 999
 999  end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getgfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
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
 501    format(A5,5(1X,1PE17.10))
      goto 10 !loop back to read statement
 11   continue !end loop when EOF is reached
c      write(0,*) "----------------------------"
             
      return
      end
