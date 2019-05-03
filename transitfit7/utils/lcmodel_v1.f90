subroutine lcmodel(nbodies,npt,tol,sol,time,itime,ans)
use precision
implicit none
integer :: nbodies,npt,nintg,i,neq,itol,itask,iopt,jt,lrw,liw,istate,j,k
integer, allocatable, dimension(:) :: iwork
real(double) :: tol,dnintg,tdnintg,t,atol,rtol,tout,tout2,jm1,tmodel,   &
 epoch
real(double), dimension(:) :: sol,time,itime,ans
real(double), allocatable, dimension(:) :: y,rwork,tcalc
real(double), allocatable, dimension(:,:) :: xpos,ypos,zpos
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

!oversampling
nintg=51!11
dnintg=dble(nintg)
tdnintg=2.0d0*dnintg

!allocate array to contain masses and position,velocity for each body.
allocate(y(nbodies*6)) !need an extra space to avoid segs (?!)
do i=1,nbodies*6
   y(i)=0.0d0
enddo
epoch=time(1)
call aei2xyz(nbodies,sol,y,epoch)
!read(5,*)

!arrays to contaim time stamps and x,y,z positions of the bodies
allocate(xpos(nbodies,nintg),ypos(nbodies,nintg),zpos(nbodies,nintg), &
   tcalc(nintg))

!setting for LSODA
lrw=22+(nbodies*6)*max(16,(nbodies*6)+9)
liw=20+6*nbodies
allocate(iwork(liw))
allocate(rwork(lrw))

nbod=nbodies
neq=6*nbodies !number of equations 6 times number of particles
t=time(1)/TU*86400.0d0     !time zero
itol=1        !atol is a scalar - for LSODA
rtol=tol*1.0d-0     !relative tolerance
atol=tol*1.0d-0     !absolute tolerance
itask=1       !don't care
iopt=0        !no advanced options
jt=2          !no jacobian, so don't worry about this one.

do i=1,npt

   tout=(time(i)-itime(i)*(0.5d0-1.0d0/tdnintg))/TU*86400.0d0

!  Deal with big gaps in the data
!   do while(tout-t.gt.itime(i)*1.5d0/TU*86400.0d0)
!         tout2=t+itime(i)/dnintg/TU*86400.0d0
!         istate=1 !error check initialization
!         call dlsoda(f, neq, y, t, tout2, itol, rtol, atol, itask,  &
!            istate, iopt, rwork, lrw, iwork, liw, jac, jt)
!         if(istate.le.0)  write(0,*) "istate: ",istate
!!         write(0,*) i,t,tout
!   enddo

!   write(0,*) 'tout:',t,tout
!   write(6,501) t,(y(6*k-5),y(6*k-4),y(6*k-3),k=1,nbodies)
   istate=1 !error check initialization
   call dlsoda(f, neq, y, t, tout, itol, rtol, atol, itask,  &
      istate, iopt, rwork, lrw, iwork, liw, jac, jt)
   if(istate.le.0)  write(0,*) "istate: ",istate
!   write(6,501) tout,(y(6*k-5),y(6*k-4),y(6*k-3),k=1,nbodies)
!   read(5,*)

   tcalc(1)=tout*TU/86400.0d0
   do j=1,nbodies
      xpos(j,1)=y(6*j-5) !X
      ypos(j,1)=y(6*j-4) !Y
      zpos(j,1)=y(6*j-3) !Z
   enddo
!   t=tout
   do j=2,nintg
      jm1=dble(j-1)
      istate=1
      tout=(time(i)-itime(i)*(0.5d0-1.0d0/tdnintg-jm1/dnintg))/TU*86400.0d0
      istate=1 !error check initialization
      call dlsoda(f, neq, y, t, tout, itol, rtol, atol, itask,  &
         istate, iopt, rwork, lrw, iwork, liw, jac, jt)
      if(istate.le.0)  write(0,*) "istate: ",istate
!      t=tout
      tcalc(j)=tout*TU/86400.0d0
      do k=1,nbodies
         xpos(k,j)=y(6*k-5) !X
         ypos(k,j)=y(6*k-4) !Y
         zpos(k,j)=y(6*k-3) !Z
      enddo
   enddo
   write(6,501) time(i),(xpos(k,2),ypos(k,2),zpos(k,2),k=2,nbodies)
   501  format(28(1PE13.6,1X))

   call transitmodel(nbodies,nintg,xpos,ypos,zpos,sol,tmodel)
   ans(i)=tmodel

enddo

return
end
