CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getfitpars(nunit,nfit,nplanet,sol,serr,err)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nunit,i,nplanet,np,j,nn,ns
      double precision sol(nfit),serr(nfit,2),err(nfit)
      character*160 command
      
c      write(0,*) "----------------------------"
c      write(0,*) "Reading in starting solution"
c      write(0,*) "----------------------------"
      
      nplanet=0 !initialize number of planets
      nn=11 !number of parameters/planet 
      ns=15 !number of parameters for the star
      
      i=0 !initialize line counter
C     Start of loop to read in file
  10  read(nunit,500,end=11,err=901) command
c        write(0,*) command
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
        elseif(command(1:3).eq."TFF") then
            read(command(5:160),*) sol(9),serr(9,1),serr(9,2),err(9)
        elseif(command(1:3).eq."LGG") then
            read(command(5:160),*) sol(10),serr(10,1),serr(10,2),err(10)
        elseif(command(1:3).eq."FEH") then
            read(command(5:160),*) sol(11),serr(11,1),serr(11,2),err(11)
        elseif(command(1:3).eq."FFF") then
            read(command(5:160),*) sol(12),serr(12,1),serr(12,2),err(12)
        elseif(command(1:3).eq."OBL") then
            read(command(5:160),*) sol(13),serr(13,1),serr(13,2),err(13)
        elseif(command(1:3).eq."OMG") then
            read(command(5:160),*) sol(14),serr(14,1),serr(14,2),err(14)
        elseif(command(1:3).eq."BET") then
            read(command(5:160),*) sol(15),serr(15,1),serr(15,2),err(15)
        elseif(command(1:2).eq."EP") then
            read(command(3:3),*,err=901) np
            if(np.gt.nplanet)nplanet=np
C           j=planet parameters*(np-1)+8 initial parameters
            j=nn*(np-1)+ns+1
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."PE") then
            read(command(3:3),*,err=901) np
            j=nn*(np-1)+ns+2
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."BB") then
            read(command(3:3),*,err=901) np
            j=nn*(np-1)+ns+3
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."RD") then
            read(command(3:3),*,err=901) np
            j=nn*(np-1)+ns+4
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."EC") then
            read(command(3:3),*,err=901) np
            j=nn*(np-1)+ns+5
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."ES") then
            read(command(3:3),*,err=901) np
            j=nn*(np-1)+ns+6
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."KR") then
            read(command(3:3),*,err=901) np
            j=nn*(np-1)+ns+7
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."TE") then
            read(command(3:3),*,err=901) np
            j=nn*(np-1)+ns+8
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."EL") then
            read(command(3:3),*,err=901) np
            j=nn*(np-1)+ns+9
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."AL") then
            read(command(3:3),*,err=901) np
            j=nn*(np-1)+ns+10
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."RT") then
            read(command(3:3),*,err=901) np
            j=nn*(np-1)+ns+11
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