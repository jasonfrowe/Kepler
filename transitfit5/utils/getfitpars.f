CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine defaultpars(nfit,nplanet,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nplanet
      real*8 sol(nfit)
      
      nplanet=2
      
      sol(1)=1.0 !mean stellar density (g/cm^3)
      sol(2)=0.0               !limb-darkening
      sol(3)=1.4638526569E-01
      sol(4)=0.0
      sol(5)=1.7078213371E-01
      sol(6)=0.0 !dilution parameter

      sol(7)=67.0 !T0 (days)
      sol(8)=3.0 !period (days)
      sol(9)=0.5 !impact parameter
      sol(10)=0.01 !Rp/R*
      sol(11)=0.0 !zero point
      sol(12)=0.0 !ecosw
      sol(13)=0.0 !esinw
      sol(14)=10.0 !RV amplitude (m/s)
      sol(15)=0.0 !RV offset
      
      sol(16)=67.0 !T0 (days)
      sol(17)=10.0 !period (days)
      sol(18)=0.5 !impact parameter
      sol(19)=0.01 !Rp/R*
      sol(20)=0.0 !zero point
      sol(21)=0.0 !ecosw
      sol(22)=0.0 !esinw
      sol(23)=10.0 !RV amplitude (m/s)
      sol(24)=0.0 !RV offset     
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getfitpars(nunit,nfit,nplanet,sol,serr,err)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nunit,i,nplanet,np,j
      double precision sol(nfit),serr(nfit,2),err(nfit)
      character*160 command

c      write(0,*) "----------------------------"
c      write(0,*) "Reading in starting solution"
c      write(0,*) "----------------------------"

      nplanet=0 !initialize number of planets

      i=0 !initialize line counter
C     Start of loop to read in file
  10  read(nunit,500,end=11,err=901) command
 500  format(A160) !command much be contained within first 80 characters
        i=i+1 !increase line counter
        if(command(1:1).eq."#") then
c            write(0,*) "Comment line on line ",i
        elseif(command(1:3).eq."RHO") then
            read(command(5:160),*) sol(1),serr(1,1),serr(1,2),err(1)
c            write(0,*) "WTF!",sol(1),serr(1,1)
        elseif(command(1:3).eq."NL1") then
            read(command(5:160),*) sol(2),serr(2,1),serr(2,2),err(2)
        elseif(command(1:3).eq."NL2") then
            read(command(5:160),*) sol(3),serr(3,1),serr(3,2),err(3)
        elseif(command(1:3).eq."NL3") then
            read(command(5:160),*) sol(4),serr(4,1),serr(4,2),err(4)
        elseif(command(1:3).eq."NL4") then
            read(command(5:160),*) sol(5),serr(5,1),serr(5,2),err(5)
        elseif(command(1:3).eq."DIL") then
            read(command(5:160),*) sol(6),serr(6,1),serr(6,2),err(6)
        elseif(command(1:3).eq."VOF") then
            read(command(5:160),*) sol(7),serr(7,1),serr(7,2),err(7)
        elseif(command(1:3).eq."ZPT") then
            read(command(5:160),*) sol(8),serr(8,1),serr(8,2),err(8)
        elseif(command(1:2).eq."EP") then
            read(command(3:3),*,err=901) np
            if(np.gt.nplanet)nplanet=np
C           j=planet parameters*(np-1)+8 initial parameters
            j=10*(np-1)+8+1
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."PE") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+2
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."BB") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+3
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."RD") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+4
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."EC") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+5
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."ES") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+6
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."KR") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+7
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."TE") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+8
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."EL") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+9
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."AL") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+10
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        endif
 501    format(A5,5(1X,1PE17.10))
c        write(0,*) command
      goto 10 !loop back to read statement
 11   continue !end loop when EOF is reached
c      write(0,*) "----------------------------"

      goto 999
 901  write(0,*) "Error on line ",i+1
      pause
      goto 999
 999  return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getfitpars2(nunit,nfit,nplanet,sol,serr,err)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nunit,i,nplanet,np,j
      double precision sol(nfit),serr(nfit,2),err(nfit,2)
      character*160 command
      
c      write(0,*) "----------------------------"
c      write(0,*) "Reading in starting solution"
c      write(0,*) "----------------------------"
      
      nplanet=0 !initialize number of planets
      
      i=0 !initialize line counter
C     Start of loop to read in file
  10  read(nunit,500,end=11,err=901) command
 500  format(A160) !command much be contained within first 80 characters
        i=i+1 !increase line counter
        if(command(1:1).eq."#") then
c            write(0,*) "Comment line on line ",i
        elseif(command(1:3).eq."RHO") then
            read(command(5:160),*) sol(1),serr(1,1),serr(1,2),err(1,1),
     .         err(1,2)
c            write(0,*) "WTF!",sol(1),serr(1,1)
        elseif(command(1:3).eq."NL1") then
            read(command(5:160),*) sol(2),serr(2,1),serr(2,2),err(2,1),
     .         err(2,2)
        elseif(command(1:3).eq."NL2") then
            read(command(5:160),*) sol(3),serr(3,1),serr(3,2),err(3,1),
     .         err(3,2)
        elseif(command(1:3).eq."NL3") then
            read(command(5:160),*) sol(4),serr(4,1),serr(4,2),err(4,1),
     .         err(4,2)
        elseif(command(1:3).eq."NL4") then
            read(command(5:160),*) sol(5),serr(5,1),serr(5,2),err(5,1),
     .         err(5,2)
        elseif(command(1:3).eq."DIL") then
            read(command(5:160),*) sol(6),serr(6,1),serr(6,2),err(6,1),
     .         err(6,2)
        elseif(command(1:3).eq."VOF") then
            read(command(5:160),*) sol(7),serr(7,1),serr(7,2),err(7,1),
     .         err(7,2)
        elseif(command(1:3).eq."ZPT") then
            read(command(5:160),*) sol(8),serr(8,1),serr(8,2),err(8,1),
     .         err(8,2)
        elseif(command(1:2).eq."EP") then
            read(command(3:3),*,err=901) np
            if(np.gt.nplanet)nplanet=np
C           j=planet parameters*(np-1)+8 initial parameters
            j=10*(np-1)+8+1
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j,1),
     .         err(j,2)
        elseif(command(1:2).eq."PE") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+2
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j,1),
     .         err(j,2)
        elseif(command(1:2).eq."BB") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+3
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j,1),
     .         err(j,2)
        elseif(command(1:2).eq."RD") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+4
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j,1),
     .         err(j,2)
        elseif(command(1:2).eq."EC") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+5
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j,1),
     .         err(j,2)
        elseif(command(1:2).eq."ES") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+6
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j,1),
     .         err(j,2)
        elseif(command(1:2).eq."KR") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+7
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j,1),
     .         err(j,2)
        elseif(command(1:2).eq."TE") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+8
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j,1),
     .         err(j,2)
        elseif(command(1:2).eq."EL") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+9
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j,1),
     .         err(j,2)
        elseif(command(1:2).eq."AL") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+10
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j,1),
     .         err(j,2)
        endif
 501    format(A5,5(1X,1PE17.10))
c        write(0,*) command
      goto 10 !loop back to read statement
 11   continue !end loop when EOF is reached
c      write(0,*) "----------------------------"

      goto 999
 901  write(0,*) "Error on line ",i+1
      pause
      goto 999       
 999  return
      end
