subroutine lcmodel(nbodies,npt,tol,sol,time,itime,ans)
use precision
implicit none
integer :: nbodies,npt,nintg,i,neq,itol,itask,iopt,jt,lrw,liw,istate,j,k
integer, allocatable, dimension(:) :: iwork,tc
real(double) :: tol,dnintg,tdnintg,t,atol,rtol,tout,tout2,jm1,tmodel,   &
 epoch
real(double), dimension(:) :: sol,time,itime,ans
real(double), allocatable, dimension(:) :: y,rwork,tcalc,opos
real(double), allocatable, dimension(:,:) :: xpos,ypos,zpos,told,bold
!mercury vars
integer :: nclo
integer, parameter :: cmax=50
integer, dimension(cmax) :: iclo,jclo
real(double), dimension(cmax) :: dclo,tclo
real(double), dimension(6,cmax) :: ixvclo,jxvclo
integer :: algor,nbig,ngflag,opflag,colflag,dtflag
integer, dimension(8) :: opt
integer, allocatable, dimension(:) :: stat
real(double) :: h0,rcen,Pi,rmax,cefac,tstart
real(double), dimension(3) :: jcen,en,am
real(double), allocatable, dimension(:) :: rho,rceh,rphys,rce,rcrit
real(double), allocatable, dimension(:,:) :: xh,vh,s,x,v
external f,jac

interface
   subroutine aei2xyz(nbodies,sol,y,epoch)
      use precision
      implicit none
      integer, intent(inout) :: nbodies
      real(double), dimension(:), intent(inout) :: sol
      real(double), dimension(:), intent(inout) :: y
      real(double), intent(inout) :: epoch
   end subroutine aei2xyz
end interface

interface
   subroutine transitmodel(nbodies,nintg,xpos,ypos,zpos,sol,tmodel)
      use precision
      implicit none
      integer, intent(in) :: nbodies,nintg
      real(double), dimension(:,:), intent(in) :: xpos,ypos,zpos
      real(double), dimension(:), intent(in) :: sol
      real(double), intent(out) :: tmodel
   end subroutine transitmodel
end interface

interface
   subroutine octiming(nbodies,nintg,xpos,ypos,zpos,sol,opos,tc,tcalc,  &
    told,bold)
      use precision
      implicit none
      integer, intent(in) :: nbodies,nintg
      integer, dimension(:), intent(inout) :: tc
      real(double), dimension(:,:), intent(in) :: xpos,ypos,zpos
      real(double), dimension(:), intent(in) :: sol,tcalc
      real(double), dimension(:), intent(inout) :: opos
      real(double), dimension(:,:), intent(inout) :: told,bold
   end subroutine octiming
end interface

!oversampling
nintg=11
dnintg=dble(nintg)
tdnintg=2.0d0*dnintg

allocate(opos(nbodies),tc(nbodies)) !used by octiming for previous position
allocate(tcalc(nintg)) !storing model times for octiming
allocate(told(nbodies,5),bold(nbodies,5))
opos=0.0d0 !initilization to zero

!allocate array to contain masses and position,velocity for each body.
allocate(y(nbodies*6)) !need an extra space to avoid segs (?!)
do i=1,nbodies*6
   y(i)=0.0d0
enddo
epoch=time(1)
call aei2xyz(nbodies,sol,y,epoch)
!read(5,*)

!arrays to contaim time stamps and x,y,z positions of the bodies
allocate(xpos(nbodies,nintg),ypos(nbodies,nintg),zpos(nbodies,nintg))

!setting for LSODA
!lrw=22+(nbodies*6)*max(16,(nbodies*6)+9)
!liw=20+6*nbodies
!allocate(iwork(liw))
!allocate(rwork(lrw))

nbod=nbodies
neq=6*nbodies !number of equations 6 times number of particles
t=time(1)!/365.25!/TU*86400.0d0     !time zero
itol=1        !atol is a scalar - for LSODA
rtol=tol*1.0d-0     !relative tolerance
atol=tol*1.0d-0     !absolute tolerance
itask=1       !don't care
iopt=0        !no advanced options
jt=2          !no jacobian, so don't worry about this one.

!Mercury vars
algor=10 !algorithm being used.
jcen=0.0d0 !J2,J4,J6 of central body
Pi=acos(-1.d0)!define Pi
!radius of star
rcen=(sol(12)/(4.0d0/3.0d0*Pi*sol(1)*1000.0d0)*Mearth)**(1.0d0/3.0d0)/AU
!write(0,*) "rcen: ",rcen
rmax=100.0d0 !distance at which particles are considered ejected (AU)
cefac=3.0d0 !Hill factor radius for close encounters (I think?)
nbig=nbod !number of 'big' bodies
allocate(xh(3,nbod),vh(3,nbod),s(3,nbod),rho(nbod),rceh(nbod),          &
 rphys(nbod),rce(nbod),rcrit(nbod),x(3,nbod),v(3,nbod),stat(nbod))
s=0.0d0 !spin angular momentum
rho=5.0d0 !density of objects. (g/cc)
rho(1)=sol(1) !density of central star (might as well use fit value)
do i=1,nbod
   rho(i)=rho(i)*1.4959787e13**3.0d0*2.959122082855911d-4/1.9891e33
enddo
rceh=0.0d0 !close-encounter limit (Hill radii)
ngflag=0 !non-grav forces is turned off.
opt=0 !options
en=0 !energies
am=0 !angular momentums
stat=0 !0-alive,1-to be removed
opflag=0 !integration mode
colflag=0 !collision flag? probably doesn't need to be set

do i=1,npt

   !target integration time.. (first step is negative!)
   tout=(time(i)-itime(i)*(0.5d0-1.0d0/tdnintg))!/365.25!/TU*86400.0d0

!   write(0,*) 'tout:',t,tout
!   write(6,501) t,(y(6*k-5),y(6*k-4),y(6*k-3),k=1,nbodies)
   h0=tout-t
   if(i.eq.1) then
      do j=1,nbodies
         xh(1,j)=y(6*j-5)
         xh(2,j)=y(6*j-4)
         xh(3,j)=y(6*j-3)
         vh(1,j)=y(6*j-2)
         vh(2,j)=y(6*j-1)
         vh(3,j)=y(6*j)
!         write(0,501) vh(1,j),vh(2,j),vh(3,j)
      enddo
!      read(5,*)
      call mce_init (t,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,       &
       m,xh,vh,s,rho,rceh,rphys,rce,rcrit,1)
      call mco_h2dh (t,jcen,nbod,nbig,h0,m,xh,vh,x,v)
      dtflag=0 !first time calling mdt_hy
      tstart=t
   else
      dtflag=2 !normal call.
   endif
   call mdt_hy (t,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,           &
      nbig,m,x,v,s,rphys,rcrit,rce,stat,algor,opt,dtflag,        &
      ngflag,opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo)
   t=t+h0
   call mco_dh2h (t,jcen,nbod,nbig,h0,m,x,v,xh,vh)
   do j=1,nbodies
      y(6*j-5)=xh(1,j)
      y(6*j-4)=xh(2,j)
      y(6*j-3)=xh(3,j)
      y(6*j-2)=vh(1,j)
      y(6*j-1)=vh(2,j)
      y(6*j)=vh(3,j)
   enddo

!   write(6,501) t,(y(6*k-5),y(6*k-4),y(6*k-3),k=1,nbodies)
!   read(5,*)

   tcalc(1)=t
   do j=1,nbodies
      xpos(j,1)=y(6*j-5) !X
      ypos(j,1)=y(6*j-4) !Y
      zpos(j,1)=y(6*j-3) !Z
   enddo
!   t=tout
   do j=2,nintg
      jm1=dble(j-1)
      istate=1
      tout=(time(i)-itime(i)*(0.5d0-1.0d0/tdnintg-jm1/dnintg))!/365.25!/TU*86400.0d0

      h0=tout-t
!      do k=1,nbodies
!         xh(1,k)=y(6*k-5)
!         xh(2,k)=y(6*k-4)
!         xh(3,k)=y(6*k-3)
!         vh(1,k)=y(6*k-2)
!         vh(2,k)=y(6*k-1)
!         vh(3,k)=y(6*k)
!      enddo
!      call mce_init (t,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,       &
!       m,xh,vh,s,rho,rceh,rphys,rce,rcrit,1)
!      call mco_h2dh (t,jcen,nbod,nbig,h0,m,xh,vh,x,v)
      dtflag=2 !normal call.
!      tstart=t
      call mdt_hy (t,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,           &
       nbig,m,x,v,s,rphys,rcrit,rce,stat,algor,opt,dtflag,        &
       ngflag,opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo)
      t=t+h0
      call mco_dh2h (t,jcen,nbod,nbig,h0,m,x,v,xh,vh)
      do k=2,nbodies
         y(6*k-5)=xh(1,k)
         y(6*k-4)=xh(2,k)
         y(6*k-3)=xh(3,k)
         y(6*k-2)=vh(1,k)
         y(6*k-1)=vh(2,k)
         y(6*k)=vh(3,k)
      enddo

      tcalc(j)=t
      do k=1,nbodies
         xpos(k,j)=y(6*k-5) !X
         ypos(k,j)=y(6*k-4) !Y
         zpos(k,j)=y(6*k-3) !Z
      enddo
   enddo
!   write(0,*) time(i),t,tout,h0
!   write(0,*)
!   write(6,501) time(i),(xpos(k,2),ypos(k,2),zpos(k,2),k=2,nbodies)
   501  format(28(1PE13.6,1X))
!   read(5,*)

!if you want the model TTVs, this routine will calculate them.
!   call octiming(nbodies,nintg,xpos,ypos,zpos,sol,opos,tc,tcalc,told,bold)

   call transitmodel(nbodies,nintg,xpos,ypos,zpos,sol,tmodel)
   ans(i)=tmodel

enddo

return
end
