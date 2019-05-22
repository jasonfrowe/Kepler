subroutine readinputsol(inputsol,nbodies,sol,serr)
use precision
implicit none
integer :: nbodies,nunit,filestatus,rnbody,i,j,np
real(double) :: rsol
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: serr
real(double), dimension(3) :: eread
character(80) :: inputsol,command

nunit=10
open(unit=nunit,file=inputsol,iostat=filestatus,status='old')

if(filestatus>0)then
   write(0,*) "Cannot open ",inputsol
   stop
endif

do
   read(nunit,*,iostat=filestatus) command,rsol,(eread(j),j=1,3)
   if(filestatus == 0) then
!      write(0,*) command
      if(command.eq.'RHO')then
         sol(1)=rsol
         serr(1,1)=eread(1)
         serr(1,2)=eread(2)
!         err(1)=eread(3)
      elseif(command.eq.'NL1')then
         sol(2)=rsol
         serr(2,1)=eread(1)
         serr(2,2)=eread(2)
!         err(2)=eread(3)
      elseif(command.eq.'NL2')then
         sol(3)=rsol
         serr(3,1)=eread(1)
         serr(3,2)=eread(2)
!         err(3)=eread(3)
      elseif(command.eq.'NL3')then
         sol(4)=rsol
         serr(4,1)=eread(1)
         serr(4,2)=eread(2)
!         err(4)=eread(3)
      elseif(command.eq.'NL4')then
         sol(5)=rsol
         serr(5,1)=eread(1)
         serr(5,2)=eread(2)
!         err(5)=eread(3)
      elseif(command.eq.'DIL')then
         sol(6)=rsol
         serr(6,1)=eread(1)
         serr(6,2)=eread(2)
!         err(6)=eread(3)
      elseif(command.eq.'ZPT')then
         sol(7)=rsol
         serr(7,1)=eread(1)
         serr(7,2)=eread(2)
!         err(7)=eread(3)
!         write(0,*) "ZPT: ",sol(8)
      else
         read(command(3:3),*) rnbody
         np=7+7*(rnbody-1)
         if(command(1:2).eq.'EP')then
            sol(np+1)=rsol
            serr(np+1,1)=eread(1)
            serr(np+1,2)=eread(2)
         elseif(command(1:2).eq.'PE')then
            sol(np+2)=rsol
            serr(np+2,1)=eread(1)
            serr(np+2,2)=eread(2)
         elseif(command(1:2).eq.'BB')then
            sol(np+3)=rsol
            serr(np+3,1)=eread(1)
            serr(np+3,2)=eread(2)
         elseif(command(1:2).eq.'RD')then
            sol(np+4)=rsol
            serr(np+4,1)=eread(1)
            serr(np+4,2)=eread(2)
         elseif(command(1:2).eq.'MA')then
            sol(np+5)=rsol
            serr(np+5,1)=eread(1)
            serr(np+5,2)=eread(2)
         elseif(command(1:2).eq.'EC')then !this is sqrt(e)cosw
            sol(np+6)=rsol
            serr(np+6,1)=eread(1)
            serr(np+6,2)=eread(2)
         elseif(command(1:2).eq.'ES')then !this is sqrt(e)sinw
            sol(np+7)=rsol
            serr(np+7,1)=eread(1)
            serr(np+7,2)=eread(2)
         endif
      endif
      cycle
   elseif(filestatus == -1) then
      exit
   else
      write(0,*) "File Error!!"
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo

close(nunit)

!do i=1,nbodies
!   m(i)=m(i)/MU*Mearth
!   y(6*i-5)=y(6*i-5)/LU*AU
!   y(6*i-4)=y(6*i-4)/LU*AU
!   y(6*i-3)=y(6*i-3)/LU*AU
!   y(6*i-2)=y(6*i-2)*TU/LU
!   y(6*i-1)=y(6*i-1)*TU/LU
!   y(6*i)=y(6*i)*TU/LU
!enddo

return
end subroutine
