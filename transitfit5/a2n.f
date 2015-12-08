      program a2n
      implicit none
      integer nfitin,nunit,iargc,nfit,nplanet,i,col
      parameter(nfitin=18,nfit=108)
      double precision solin(nfitin),serrin(nfitin,2),errin(nfitin),
     .  sol(nfit),serr(nfit,2),err(nfit),Pi,tPi,pid2,G,Per,Psec,eccn,
     .  adrs,rhos
      character*3 titles(nfit)
      character*80 inputsol
      
      do 10 i=1,nfit !initialize values
        serr(i,1)=0.0
        serr(i,2)=0.0
        err(i)=0.0
 10   continue 
      
      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi 
      pid2=Pi/2.0d0
      G=6.674d-11 !N m^2 kg^-2  Gravitation constamt
      
      if(iargc().lt.1) goto 901
      
      call getarg(1,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfitin,solin,serrin,errin)
      close(nunit) !release unit number as we are done with file
      
      nplanet=1
      
      per=solin(2)     !Period (days)   
      Psec=per*24.0d0*60.0d0*60.0d0 !period (sec)   
      eccn=sqrt(solin(7)*solin(7)+solin(8)*solin(8)) !eccentricity
      if(eccn.ge.1.0) eccn=0.99 
C     a/R*
      adrs=solin(5)*per/tpi*sqrt(1-solin(3))*(1+solin(8))/
     .  sqrt(1-eccn*eccn)
      rhos= adrs**3.0*Pi*3.0d0/(Psec*Psec*G) !kg/m^3
      
      sol(1)=rhos/1000.0 !g/cm^3
      serr(1,2)=-1.0
      sol(2)=solin(11) !limb-darkening
      sol(3)=solin(12)
      sol(4)=solin(13)
      sol(5)=solin(14)
      sol(6)=solin(18) !dilution
      sol(7)=solin(10) !gamma (m/s)
      sol(8)=solin(6)  !zpt
      serr(8,2)=-1.0
      sol(9)=solin(1) !EPO (BJD-2454900)
      serr(9,2)=-1.0
      sol(10)=solin(2) !Per (days)
      serr(10,2)=-1.0
      sol(11)=sqrt(solin(3)) !b^2
      serr(11,2)=-1.0
      sol(12)=solin(4) !r/R*
      serr(12,2)=-1.0
      sol(13)=solin(7) !ECW
      sol(14)=solin(8) !ESW
      sol(15)=solin(9) !K (m/s)
      sol(16)=0.0!solin(15) !TED (ppm)
      sol(17)=solin(16) !Ellip (ppm)
      sol(18)=solin(17) !albedo (ppm)
      
      do 11 i=2,iargc()
        call getarg(i,inputsol) !get filename for input solution
        nunit=10 !unit number used for file input
        open(unit=nunit,file=inputsol,status='old',err=902)
C       We start by reading in solution from input file
        call getfitpars(nunit,nfitin,solin,serrin,errin)
        close(nunit) !release unit number as we are done with file
        
        nplanet=nplanet+1
        col=8+10*(nplanet-1)
        
        sol(col+1)=solin(1) !EPO
        serr(col+1,2)=-1.0
        sol(col+2)=solin(2) !Per
        serr(col+2,2)=-1.0
        sol(col+3)=solin(3) !b^2
        serr(col+3,2)=-1.0
        sol(col+4)=solin(4) !r/R*
        serr(col+4,2)=-1.0
        sol(col+5)=solin(7) !ecw
        sol(col+6)=solin(8) !esw
        sol(col+7)=solin(9) !K (m/s)
        sol(col+8)=0.0!solin(15) !TED
        sol(col+9)=solin(16) !ELL
        sol(col+10)=solin(17) !ALB
 11   continue
      
      call exportfit(nfit,nplanet,sol,serr,err,titles)
      
      goto 999
 901  write(0,*) "Usage: a2n a1.dat a2.dat a3.dat..."
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exportfit(nfit,nplanet,sol,serr,err,titles)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,i,nplanet
      double precision sol(nfit),serr(nfit,2),err(nfit)
      character*3 tit
      include 'titles.f'
      
      open(unit=11,file="newfit.dat")
      
      do 10 i=1,nplanet*10+8
        if(i.le.8)then
            tit=titles(i)
        else
            tit=titles(i-10*((i-9)/10))
            write(tit(3:3),501) ((i-9)/10)+1
 501        format(I1)
c            write(0,*) i,i-10*((i-9)/10)
        endif
        write(11,500) tit,sol(i),serr(i,1),serr(i,2),err(i)
 10   continue
 500  format(A3,5(1X,1PE17.10))
 
      close(11)
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getfitpars(nunit,nfit,sol,serr,err)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nunit,i
      double precision sol(nfit),serr(nfit,2),err(nfit)
      character*160 command
      
c      write(0,*) "----------------------------"
c      write(0,*) "Reading in starting solution"
c      write(0,*) "----------------------------"
      
      i=0 !initialize line counter
C     Start of loop to read in file
  10  read(nunit,500,end=11) command
 500  format(A160) !command much be contained within first 80 characters
        i=i+1 !increase line counter
        if(command(1:1).eq."#") then
            write(0,*) "Comment line on line ",i
        elseif(command(1:3).eq."EPO") then
            read(command(5:160),*) sol(1),serr(1,1),serr(1,2),err(1)
        elseif(command(1:3).eq."PER") then
            read(command(5:160),*) sol(2),serr(2,1),serr(2,2),err(2)
        elseif(command(1:3).eq."BBB") then
            read(command(5:160),*) sol(3),serr(3,1),serr(3,2),err(3)
        elseif(command(1:3).eq."RDR") then
            read(command(5:160),*) sol(4),serr(4,1),serr(4,2),err(4)
        elseif(command(1:3).eq."ZDR") then
            read(command(5:160),*) sol(5),serr(5,1),serr(5,2),err(5)
        elseif(command(1:3).eq."ZPT") then
            read(command(5:160),*) sol(6),serr(6,1),serr(6,2),err(6)
        elseif(command(1:3).eq."ECW") then
            read(command(5:160),*) sol(7),serr(7,1),serr(7,2),err(7)
        elseif(command(1:3).eq."ESW") then
            read(command(5:160),*) sol(8),serr(8,1),serr(8,2),err(8)
        elseif(command(1:3).eq."KRV") then
            read(command(5:160),*) sol(9),serr(9,1),serr(9,2),err(9)
        elseif(command(1:3).eq."VOF") then
            read(command(5:160),*) sol(10),serr(10,1),serr(10,2),err(10)
        elseif(command(1:3).eq."NL1") then
            read(command(5:160),*) sol(11),serr(11,1),serr(11,2),err(11)
        elseif(command(1:3).eq."NL2") then
            read(command(5:160),*) sol(12),serr(12,1),serr(12,2),err(12)
        elseif(command(1:3).eq."NL3") then
            read(command(5:160),*) sol(13),serr(13,1),serr(13,2),err(13)
        elseif(command(1:3).eq."NL4") then
            read(command(5:160),*) sol(14),serr(14,1),serr(14,2),err(14)
        elseif(command(1:3).eq."TED") then
            read(command(5:160),*) sol(15),serr(15,1),serr(15,2),err(15)
        elseif(command(1:3).eq."ELL") then
            read(command(5:160),*) sol(16),serr(16,1),serr(16,2),err(16)
        elseif(command(1:3).eq."ALB") then
            read(command(5:160),*) sol(17),serr(17,1),serr(17,2),err(17)
        elseif(command(1:3).eq."DIL") then
            read(command(5:160),*) sol(18),serr(18,1),serr(18,2),err(18)
        endif
 501    format(A5,5(1X,1PE17.10))
      goto 10 !loop back to read statement
 11   continue !end loop when EOF is reached
c      write(0,*) "----------------------------"
             
      return
      end 
