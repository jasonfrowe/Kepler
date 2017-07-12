CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine defaultpars(nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit
      real*8 sol(nfit)
      
      sol(1)=55.90128   !center of transit time (days)
      sol(2)=3.548460d0 !Period (days)
      sol(3)=0.584      !impact parameter (b)
      sol(4)=0.0828     !Rp/R* 
      sol(5)=10.81       !z/R*
      sol(6)=0.0        !zeropoint     
      sol(7)=0.0d0      !ecosw
      sol(8)=0.0d0      !esinw
      sol(9)=0.0d0      !K (RV amplitude)
      sol(10)=0.0d0     !RV offset
      sol(11)=0.0       !limb darkening
      sol(12)=1.4638526569E-01
      sol(13)=0.0
      sol(14)=1.7078213371E-01
      sol(15)=0.0        !Depth of secondary
      sol(16)=0.0        !amplitude of ellipsodial variations
      sol(17)=0.0       !amplitude of phase changes (e.g. Albedo)
      sol(18)=0.0       !dilution parameter
      
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
 501    format(A5,5(1X,PE17.10))
      goto 10 !loop back to read statement
 11   continue !end loop when EOF is reached
c      write(0,*) "----------------------------"
             
      return
      end 