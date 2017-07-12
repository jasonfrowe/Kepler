CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine defaultpars(nfit,nplanet,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nplanet
      real*8 sol(nfit)
      
      nplanet=1
      
      sol(1)=1.0 !stellar mass
      sol(2)=1.0 !stellar radius
      sol(3)=0.1966 !limb-darkening
      sol(4)=0.2908
      sol(5)=0.0
      sol(6)=0.0
      sol(7)=0.0 !dilution
      sol(8)=0.0 !gamma
      sol(9)=0.0 !photometric zeropoint
      sol(10)=1.0 !Doppler Scaling parameter
      sol(11)=1.0 !ellipsoidal scaling
      
      sol(12)=76.069060641 !T0
      sol(13)=23.876134918 !period
      sol(14)=0.35 !b
      sol(15)=66000.0 !Mp
      sol(16)=10.0 !Rp
      sol(17)=0.0 !ecw
      sol(18)=0.0 !esw
      sol(19)=0.0 !albedo
      sol(20)=0.0 !occultation depth
      
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
        elseif(command(1:3).eq."MAS") then
            read(command(5:160),*) sol(1),serr(1,1),serr(1,2),err(1)
        elseif(command(1:3).eq."RAD") then
            read(command(5:160),*) sol(2),serr(2,1),serr(2,2),err(2)
        elseif(command(1:3).eq."NL1") then
            read(command(5:160),*) sol(3),serr(3,1),serr(3,2),err(3)
        elseif(command(1:3).eq."NL2") then
            read(command(5:160),*) sol(4),serr(4,1),serr(4,2),err(4)
        elseif(command(1:3).eq."NL3") then
            read(command(5:160),*) sol(5),serr(5,1),serr(5,2),err(5)
        elseif(command(1:3).eq."NL4") then
            read(command(5:160),*) sol(6),serr(6,1),serr(6,2),err(6)
        elseif(command(1:3).eq."DIL") then
            read(command(5:160),*) sol(7),serr(7,1),serr(7,2),err(7)
        elseif(command(1:3).eq."VOF") then
            read(command(5:160),*) sol(8),serr(8,1),serr(8,2),err(8)
        elseif(command(1:3).eq."ZPT") then
            read(command(5:160),*) sol(9),serr(9,1),serr(9,2),err(9)
        elseif(command(1:3).eq."DOP") then
            read(command(5:160),*) sol(10),serr(10,1),serr(10,2),err(10)
        elseif(command(1:3).eq."ELL") then
            read(command(5:160),*) sol(11),serr(11,1),serr(11,2),err(11)
        elseif(command(1:2).eq."EP") then
            read(command(3:3),*,err=901) np
            if(np.gt.nplanet)nplanet=np
            j=9*(np-1)+11+1
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."PE") then
            read(command(3:3),*,err=901) np
            j=9*(np-1)+11+2
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."BB") then
            read(command(3:3),*,err=901) np
            j=9*(np-1)+11+3
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."MP") then
            read(command(3:3),*,err=901) np
            j=9*(np-1)+11+4
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."RP") then
            read(command(3:3),*,err=901) np
            j=9*(np-1)+11+5
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."EC") then
            read(command(3:3),*,err=901) np
            j=9*(np-1)+11+6
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."ES") then
            read(command(3:3),*,err=901) np
            j=9*(np-1)+11+7
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."AL") then
            read(command(3:3),*,err=901) np
            j=9*(np-1)+11+8
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."TE") then
            read(command(3:3),*,err=901) np
            j=9*(np-1)+11+9
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


