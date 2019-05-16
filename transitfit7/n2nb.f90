program n2nb
use precision
implicit none
integer :: nfit,nplanet,iargc,npt,nmax,i,j
real(double) :: tPi,T0,Per,EP,Psec,incl,th1,th2,vel,M1,M2,asemi
real(double), allocatable, dimension(:) :: sol,time,flux,ferr,itime, &
   solout,err,y,yerr
real(double), allocatable, dimension(:) :: m,merr
real(double), allocatable, dimension(:,:) :: serr,yserr,mserr
character(80) :: photfile,inputsol,charmass

interface
   subroutine readkeplc(photfile,npt,time,flux,ferr,itime)
      use precision, only: double
      implicit none
      character(80), intent(in) :: photfile
      integer, intent(out) :: npt
      real(double), dimension(:), intent(inout) :: time,flux,ferr,itime
   end subroutine readkeplc
end interface

interface
   subroutine calcnbodies(inputsol,nbodies)
      use precision, only: double
      implicit none
      character(80), intent(in) :: inputsol
      integer, intent(out) :: nbodies
   end subroutine calcnbodies
end interface

interface
   subroutine getfitpars(inputsol,sol)
      use precision, only: double
      implicit none
      character(80), intent(in) :: inputsol
      real(double), dimension(:) :: sol
   end subroutine getfitpars
end interface

interface
   subroutine exportfit(nplanet,solout,serr,err,y,yerr,yserr)
      use precision, only: double
      implicit none
      integer, intent(in) :: nplanet
      real(double), dimension(:), intent(in) :: solout,err,y,yerr
      real(double), dimension(:,:), intent(in) :: serr,yserr
   end subroutine exportfit
end interface

if(iargc().lt.4) then
   write(0,*) "Usage: n2nb <photometry> <n.dat> <mstar> <mi>"
   write(0,*) " <mstar> Mass of star [Msun]"
   write(0,*) " <Mi> Mass of each planet i [Mearth]"
   stop
endif

call getarg(1,photfile)
nmax=100000
allocate(time(nmax),flux(nmax),ferr(nmax),itime(nmax))
call readkeplc(photfile,npt,time,flux,ferr,itime)
write(0,*) "Number of observations read: ",npt

call getarg(2,inputsol)
call calcnbodies(inputsol,nplanet)
write(0,*) "Number of planets", nplanet

nfit=8+nplanet*10
allocate(sol(nfit),solout(8+nplanet+1),serr(8+nplanet+1,2), &
   err(8+nplanet+1),m(nplanet+1),mserr(nplanet+1,2),merr(nplanet+1), &
   y((nplanet+1)*6),yserr((nplanet+1)*6,2),yerr((nplanet+1)*6))

!get masses
if(iargc().lt.3+nplanet) then
   write(0,*) "Usage: n2nb <photometry> <n.dat> <mstar> <mi>"
   write(0,*) " <mstar> Mass of star [Msun]"
   write(0,*) " <Mi> Mass of each planet i [Mearth]"
   write(0,*) "**Must list mass for all planets: ",nplanet
   stop
endif
call getarg(3,charmass)
read(charmass,*) m(1)
m(1)=m(1)*Msun/Mearth !convert to Earth-masses
mserr(1,1)=0.0d0
mserr(1,2)=-1.0d0
merr(1)=0.0d0
do i=1,nplanet
   call getarg(3+i,charmass)
   read(charmass,*) m(i+1)
   mserr(i+1,1)=0.0d0
   mserr(i+1,2)=-1.0d0
   merr(i+1)=0.0d0
enddo

call getfitpars(inputsol,sol)

do i=1,8
   solout(i)=sol(i)
   serr(i,1)=0.0d0
   serr(i,2)=0.0d0
   err(i)=0.0d0
enddo
serr(1,2)=-1.0d0 !fit for rhostar
serr(8,2)=-1.0d0 !fit for ZPT

solout(8+1)=0.0d0 !RDR for central star is not used.
serr(8+1,1)=0.0d0
serr(8+1,2)=0.0d0
err(8+1)=0.0d0
do i=1,nplanet
   solout(8+i+1)=sol(10*(i-1)+8+4)
   serr(8+i+1,2)=-1.0d0 !fit for R/R*
   serr(8+i+1,1)=0.0d0
   err(8+i+1)=0.0d0
!   write(6,*) 8+i+1,solout(8+i+1),sol(10*(i-1)+8+4)
enddo

!central star is stationary
do i=1,6
   y(i)=0.0d0
   yserr(i,1)=0.0d0
   yserr(i,2)=-1.0d0
   yerr(i)=0.0d0
enddo

tPi=2.0d0*Pi

T0=time(1)
M1=m(1)*MEarth
do i=1,nplanet
   Ep =sol(10*(i-1)+8+1)
   Per=sol(10*(i-1)+8+2)
   M2=M(i+1)*Mearth !convert to kg
   Psec=Per*24.0*60.0*60.0
   asemi=(Psec*Psec*Gr*(M1+M2)/(4.0*Pi*Pi))**(1.0/3.0)
   vel=tPi*asemi/Psec
   incl=0.0d0
   th1=(Ep-T0)/Per-int((Ep-T0)/Per)
   th1=tPi-th1*tPi
   th2=th1+Pi/2.0d0
   y(6*i+1)=asemi*cos(th1)/AU
   y(6*i+2)=asemi*sin(th1)*cos(incl)/AU
   y(6*i+3)=asemi*sin(th1)*sin(incl)/AU
   y(6*i+4)=vel*cos(th2)
   y(6*i+5)=vel*sin(th2)*cos(incl)
   y(6*i+6)=vel*sin(th2)*sin(incl)
   do j=1,6
      yserr(6*i+j,1)=0.0d0
      yserr(6*i+j,2)=-1.0d0
      yerr(6*i+j)=0.0d0
   enddo
!   write(0,500) m(i+1),(y(j),j=6*(i+1)-5,6*(i+1))
enddo
500  format(28(1PE10.3,1X))

call exportfit(nplanet,solout,serr,err,y,yerr,yserr)

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
