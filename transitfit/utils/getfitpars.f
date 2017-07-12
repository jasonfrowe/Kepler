CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine printsol(nfit,sol,serr,Dpvary,err,doe,toff)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,i
      double precision sol(nfit),serr(nfit,2),Dpvary(nfit),toff,
     .  err(nfit),doe
      include "utils/titles.f"
      
      do 10 i=1,nfit
        write(0,500) titles(i),":",sol(i),serr(i,1),serr(i,2),Dpvary(i),
     .      err(i)
 10   continue
 500  format(A3,A1,5(1X,1PE17.10))
c      write(0,501) "OFF:",toff,doe
 501  format(A4,2(1PE17.10,1X))
      write(0,*) "--------------------------------------------"

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine defaultpars(nfit,sol,serr,Dpvary,err,doe,toff)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,i
      double precision sol(nfit),serr(nfit,2),Dpvary(nfit),toff,
     .  err(nfit),doe
C     A fit that matches HD209458 (for something to start with)

      sol(1)=1.11d0
      sol(2)=0.00d0
      sol(3)=1.14d0
      sol(4)=1.36d0
      sol(5)=3.52d0
      sol(6)=86.9d0
      sol(7)=0.0d0
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
      
      Dpvary(1)=0.02
      Dpvary(2)=0.02
      Dpvary(3)=0.02
      Dpvary(4)=0.02
      Dpvary(5)=4.0d-7
      Dpvary(6)=0.02
      Dpvary(7)=0.02 !30 minute variation (if fail, try 45)
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
      do 10 i=1,nfit
        serr(i,1)= 0.0d0
        serr(i,2)= 0.0d0
        err(i)=0.0
 10   continue
      serr(4,2)=-1.0
      serr(5,2)=-1.0
      serr(6,2)=-1.0
      serr(7,2)=-1.0
      serr(8,2)=-1.0
      toff=0.0d0
      doe=0.0
 
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
c        	write(0,501) "SMA: ",sol(1),serr(1,1),serr(1,2),Dpvary(1),
c     .          err(1)
        elseif(command(1:3).eq."PMA") then
        	read(command(5:160),*) sol(2),serr(2,1),serr(2,2),Dpvary(2),
     .          err(2)
c        	write(0,501) "PMA: ",sol(2),serr(2,1),serr(2,2),Dpvary(2),
c     .          err(2)
        elseif(command(1:3).eq."SRA") then
        	read(command(5:160),*) sol(3),serr(3,1),serr(3,2),Dpvary(3),
     .          err(3)
c        	write(0,501) "SRA: ",sol(3),serr(3,1),serr(3,2),Dpvary(3),
c     .          err(3)
        elseif(command(1:3).eq."PRA") then
        	read(command(5:160),*) sol(4),serr(4,1),serr(4,2),Dpvary(4),
     .          err(4)
c        	write(0,501) "PRA: ",sol(4),serr(4,1),serr(4,2),Dpvary(4),
c     .          err(4)
        elseif(command(1:3).eq."PER") then
        	read(command(5:160),*) sol(5),serr(5,1),serr(5,2),Dpvary(5),
     .          err(5)
c        	write(0,501) "PER: ",sol(5),serr(5,1),serr(5,2),Dpvary(5),
c     .          err(5)
        elseif(command(1:3).eq."INC") then
        	read(command(5:160),*) sol(6),serr(6,1),serr(6,2),Dpvary(6),
     .          err(6)
c        	write(0,501) "INC: ",sol(6),serr(6,1),serr(6,2),Dpvary(6),
c     .          err(6)
        elseif(command(1:3).eq."EPO") then
        	read(command(5:160),*) sol(7),serr(7,1),serr(7,2),Dpvary(7),
     .          err(7)
c        	write(0,501) "EPO: ",sol(7),serr(7,1),serr(7,2),Dpvary(7),
c     .          err(7)
        elseif(command(1:3).eq."ZPT") then
        	read(command(5:160),*) sol(8),serr(8,1),serr(8,2),Dpvary(8),
     .          err(8)
c        	write(0,501) "ZPT: ",sol(8),serr(8,1),serr(8,2),Dpvary(8),
c     .          err(8)
        elseif(command(1:3).eq."ALB") then
        	read(command(5:160),*) sol(9),serr(9,1),serr(9,2),Dpvary(9),
     .          err(9)
c        	write(0,501) "ALB: ",sol(9),serr(9,1),serr(9,2),Dpvary(9),
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
