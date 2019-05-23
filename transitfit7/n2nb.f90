program n2nb
use precision
implicit none
integer :: nfit,nplanet,iargc,npt,nmax,i,j,npout,np,nbodies
real(double) :: c1,c2,c3,c4,q1,q2
real(double), allocatable, dimension(:) :: sol,solout
real(double), allocatable, dimension(:) :: m,merr
real(double), allocatable, dimension(:,:) :: serr,yserr,mserr
character(80) :: photfile,inputsol,charmass

interface
   subroutine calcnbodies(inputsol,nbodies)
      use precision, only: double
      implicit none
      character(80), intent(in) :: inputsol
      integer, intent(out) :: nbodies
   end subroutine calcnbodies
   subroutine getfitpars(inputsol,sol)
      use precision, only: double
      implicit none
      character(80), intent(in) :: inputsol
      real(double), dimension(:) :: sol
   end subroutine getfitpars
   subroutine exportfit(nbodies,sol,serr)
      use precision
      implicit none
      integer nbodies
      real(double), dimension(:), intent(inout) :: sol
      real(double), dimension(:,:), intent(inout) :: serr
   end subroutine exportfit
end interface

if(iargc().lt.3) then
   write(0,*) "Usage: n2nb <n.dat> <mstar> <mi>"
   write(0,*) " <mstar> Mass of star [Msun]"
   write(0,*) " <Mi> Mass of each planet i [Mearth]"
   stop
endif


call getarg(1,inputsol)
call calcnbodies(inputsol,nplanet)
nbodies=nplanet+1
write(0,*) "Number of planets", nplanet

nfit=8+nplanet*10
allocate(sol(nfit),m(nplanet+1))

!get masses
if(iargc().lt.2+nplanet) then
   write(0,*) "Usage: n2nb <photometry> <n.dat> <mstar> <mi>"
   write(0,*) " <mstar> Mass of star [Msun]"
   write(0,*) " <Mi> Mass of each planet i [Mearth]"
   write(0,*) "**Must list mass for all planets: ",nplanet
   stop
endif
call getarg(2,charmass)
read(charmass,*) m(1)
m(1)=m(1)*Msun/Mearth !convert to Earth-masses
do i=1,nplanet
   call getarg(2+i,charmass)
   read(charmass,*) m(i+1)
enddo

call getfitpars(inputsol,sol)

allocate(solout(7+nbodies*7),serr(7+nbodies*7,2)) !use nbodies to allocate
serr=0.0d0 !set default to zero.
solout=0.0d0

!mean stellar density
solout(1)=sol(1)
serr(1,2)=-1.0d0

!limb-darkening, dilution

c1=sol(2)
c2=sol(3)
c3=sol(4)
c4=sol(5)

if((c3.eq.0.0).and.(c4.eq.0.0))then
   q1=(c1+c2)*(c1+c2)
   q2=c1/(2*(c1+c2))
   solout(2)=0.0
   solout(3)=0.0
   solout(4)=q1
   solout(5)=q2
elseif((c1.eq.0.0).and.(c2.eq.0.0))then
   solout(2)=0.0
   solout(3)=0.0
   solout(4)=c3
   solout(5)=c4
else
   do i=2,6
      solout(i)=sol(i)
   enddo
endif

!photometric zero point
solout(7)=sol(8)
serr(7,2)=-1.0d0

!star
!do i=8,14
!   solout(i)=0.0
!   serr(i,1)=0.0
!   serr(i,2)=0.0
!enddo
solout(12)=m(1) !mass of star

!individual planets
do i=2,nbodies
   npout=7+7*(i-1) !nb.dat
   np=10*(i-2)+8  !n0.dat
   solout(npout+1)=sol(np+1)  !epoch
   solout(npout+2)=sol(np+2)  !period
   solout(npout+3)=sol(np+3)  !bb
   solout(npout+4)=sol(np+4)  !r/r*
   solout(npout+5)=m(i)       !mass
   if(sol(np+5).ne.0.0d0)then
      solout(npout+6)=sol(np+5)  !ecosw
   else
      solout(npout+6)=1.0d-2
   endif
   if(sol(np+6).ne.0.0d0)then
      solout(npout+7)=sol(np+6)  !esinw
   else
      solout(npout+7)=1.0d-2
   endif
   do j=1,7
      serr(npout+j,2)=-1.0d0
   enddo
enddo

500  format(28(1PE10.3,1X))

call exportfit(nbodies,solout,serr)

end program n2nb

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getfitpars(inputsol,sol)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision, only: double
implicit none
integer nunit,filestatus,j
real(double) :: rsol,rnbody
real(double), dimension(:) :: sol
character(3) :: command
character(80) :: inputsol

nunit=10
open(unit=nunit,file=inputsol,iostat=filestatus,status='old')

if(filestatus>0)then
   write(0,*) "Cannot open ",inputsol
   stop
endif

do
   read(nunit,*,iostat=filestatus) command,rsol
   if(filestatus == 0) then
      if(command.eq.'RHO')then
         sol(1)=rsol
      elseif(command.eq.'NL1')then
         sol(2)=rsol
      elseif(command.eq.'NL2')then
         sol(3)=rsol
      elseif(command.eq.'NL3')then
         sol(4)=rsol
      elseif(command.eq.'NL4')then
         sol(5)=rsol
      elseif(command.eq.'DIL')then
         sol(6)=rsol
      elseif(command.eq.'VOF')then
         sol(7)=rsol
      elseif(command.eq.'ZPT')then
         sol(8)=rsol
      else
         read(command(3:3),*) rnbody
         j=10*(rnbody-1)+8
         if(command(1:2).eq.'EP')then
            sol(j+1)=rsol
         elseif(command(1:2).eq.'PE')then
            sol(j+2)=rsol
         elseif(command(1:2).eq.'BB')then
            sol(j+3)=rsol
         elseif(command(1:2).eq.'RD')then
            sol(j+4)=rsol
         elseif(command(1:2).eq.'EC')then
            sol(j+5)=rsol
         elseif(command(1:2).eq.'ES')then
            sol(j+6)=rsol
         elseif(command(1:2).eq.'KR')then
            sol(j+7)=rsol
         elseif(command(1:2).eq.'TE')then
            sol(j+8)=rsol
         elseif(command(1:2).eq.'EL')then
            sol(j+9)=rsol
         elseif(command(1:2).eq.'AL')then
            sol(j+10)=rsol
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

return
end subroutine getfitpars
