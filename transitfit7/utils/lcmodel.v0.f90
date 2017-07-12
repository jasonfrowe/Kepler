subroutine lcmodel(nbodies,npt,tol,y,sol,time,itime,ans)
use precision
implicit none
integer :: lrw,liw,neq,itol,nbodies,i,j,k,npt,itask,istate,iopt, &
   jt,nintg
integer, allocatable, dimension(:) :: iwork
real(double), allocatable, dimension(:) :: rwork,tcalc
real(double), allocatable, dimension(:,:) :: xpos,ypos,zpos
real(double) t,rtol,atol,tol,dt,tout,dnintg,tdnintg,jm1,tmodel,told, &
   tout2
real(double), dimension(:) :: time,itime
real(double), dimension(:) :: y,ans
real(double), dimension(:) :: sol
external f,jac

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
nintg=11
dnintg=dble(nintg)
tdnintg=2.0d0*dnintg

!do i=1,nbodies
!   write(0,*) i,m(i)
!enddo
!read(5,*)
!write(0,*) "tol:",tol

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
itol=1        !tolerence time - for LSODA
rtol=tol      !relative tolerance
atol=tol      !absolute tolerance
itask=1       !don't care
iopt=0        !no advanced options
jt=2          !no jacobian, so don't worry about this one.

do i=1,npt
   istate=1 !error check initialization

   tout=(time(i)-itime(i)*(0.5d0-1.0d0/tdnintg))/TU*86400.0d0

!  Deal with big gaps in the data
   do while(tout-t.gt.itime(i)*1.5d0/TU*86400.0d0)
         tout2=t+itime(i)*1.5d0/TU*86400.0d0
         call dlsoda(f, neq, y, t, tout2, itol, rtol, atol, itask,  &
            istate, iopt, rwork, lrw, iwork, liw, jac, jt)
         if(istate.le.0)  write(0,*) "istate: ",istate
   enddo

!   write(6,501) t,tout
   call dlsoda(f, neq, y, t, tout, itol, rtol, atol, itask,  &
      istate, iopt, rwork, lrw, iwork, liw, jac, jt)
   if(istate.le.0)  write(0,*) "istate: ",istate
   tcalc(1)=tout*TU/86400.0d0
   do j=1,nbodies
      xpos(j,1)=y(6*j-5) !X
      ypos(j,1)=y(6*j-4) !Y
      zpos(j,1)=y(6*j-3) !Z
   enddo
   t=tout
   do j=2,nintg
      jm1=dble(j-1)
      istate=1
      tout=(time(i)-itime(i)*(0.5d0-1.0d0/tdnintg-jm1/dnintg))/TU*86400.0d0
      call dlsoda(f, neq, y, t, tout, itol, rtol, atol, itask,  &
         istate, iopt, rwork, lrw, iwork, liw, jac, jt)
      if(istate.le.0)  write(0,*) "istate: ",istate
      t=tout
      tcalc(j)=tout*TU/86400.0d0
      do k=1,nbodies
         xpos(k,j)=y(6*k-5) !X
         ypos(k,j)=y(6*k-4) !X
         zpos(k,j)=y(6*k-3) !Z
      enddo
   enddo
   501  format(28(1PE13.6,1X))

   call transitmodel(nbodies,nintg,xpos,ypos,zpos,sol,tmodel)

!   write(6,501) time(i),tmodel,(xpos(j,5),ypos(j,5),zpos(j,5),j=1,nbodies)

   ans(i)=tmodel

enddo

return
end subroutine
