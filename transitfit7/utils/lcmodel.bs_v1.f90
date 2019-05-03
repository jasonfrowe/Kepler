subroutine lcmodel(nbodies,npt,tol,sol,time,itime,percor,ans)
use precision
implicit none
!import vars
integer, target :: nbodies
integer :: npt
real(double) :: tol
real(double), dimension(:) :: sol,time,itime,ans,percor
!local vars
integer :: nintg,i,j,k,nb,neq
real(double) :: dnintg,tdnintg,epoch,Rmin,t,tout,h0,hdid,hnext,jm1,tmodel,eps,tsec
real(double), allocatable, dimension(:), target :: m
real(double), allocatable, dimension(:) :: y,ydot,yscal
real(double), allocatable, dimension(:,:) :: xpos,ypos,zpos
!octiming vars
integer, allocatable, dimension(:) :: tc
real(double), allocatable, dimension(:) :: opos,tcalc
real(double), allocatable, dimension(:,:) :: told,bold
!plotting stuff
integer :: iplot
real :: bb(4),rasemi,rx,ry
external f

interface
   subroutine aei2xyz(nbodies,sol,y,m,epoch,percor)
      use precision
      implicit none
      integer, intent(inout) :: nbodies
      real(double), dimension(:), intent(inout) :: sol,percor
      real(double), dimension(:), intent(inout) :: y,m
      real(double), intent(inout) :: epoch
   end subroutine aei2xyz
   subroutine transitmodel(nbodies,nintg,xpos,ypos,zpos,sol,tmodel)
      use precision
      implicit none
      integer, intent(in) :: nbodies,nintg
      real(double), dimension(:,:), intent(in) :: xpos,ypos,zpos
      real(double), dimension(:), intent(in) :: sol
      real(double), intent(out) :: tmodel
   end subroutine transitmodel
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

!plotting
iplot=0 !set iplot=1 for plotting 

!definition of small number
eps=EPSILON(1.0d0)
!eps=1.0d-6

!oversampling
nintg=11
dnintg=dble(nintg)
tdnintg=2.0d0*dnintg

!arrays for octiming 
allocate(opos(nbodies),tc(nbodies)) !used by octiming for previous position
allocate(tcalc(nintg)) !storing model times for octiming
allocate(told(nbodies,5),bold(nbodies,5))
opos=0.0d0 !initilization to zero
tc=0 !initializtion of transit conidition flag to zero.

!allocate array to contain masses and position,velocity for each body.
allocate(y(nbodies*6),m(nbodies),ydot(nbodies*6),yscal(nbodies*6))
do i=1,nbodies*6
   y(i)=0.0d0
enddo
epoch=time(1)
call aei2xyz(nbodies,sol,y,m,epoch,percor)
!read(5,*)

!plotting stuff...
if(iplot.eq.1)then 
   call pgopen('/xserve') !open PGPlot device
   call pgpap(8.0,1.0) !page size
   call pgslw(3) !thicker lines
   rasemi=0.0
   do k=1,nbodies
      rasemi=max(rasemi,real(y(6*k-5))) !X
      rasemi=max(rasemi,real(y(6*k-4))) !Y
      rasemi=max(rasemi,real(y(6*k-3))) !Z
   enddo
   bb(1)=-3.0*rasemi
   bb(2)= 3.0*rasemi
   bb(3)=-3.0*rasemi
   bb(4)= 3.0*rasemi
   call pgwindow(bb(1),bb(2),bb(3),bb(4))
   CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
   call pgslw(1) !thicker lines
endif

!arrays to contaim time stamps and x,y,z positions of the bodies
allocate(xpos(nbodies,nintg),ypos(nbodies,nintg),zpos(nbodies,nintg))

!pointers for BSSTEP function "f"
nbodies2 => nbodies
m2 => m

neq=nbodies*6 !number of coupled DEs
t=time(1)    !time zero [days]
!write(0,*) "t0: ",t
do i=1,npt
   !write(0,*) "tmid: ",time(i)
   do j=1,nintg
      !target integration time.. (first step is negative!)
      jm1=dble(j-1)
      tout=(time(i)-itime(i)*(0.5d0-1.0d0/tdnintg-jm1/dnintg)) ![days]
      !write(0,*) 'tt',time(i),tout

      h0=day*(tout-t)/TU !required time-step [scaled]

      !call to the integrator
      !write(0,*) 'h0',h0

      hdid=0.0d0
      do while(abs(h0-hdid).gt.2.0*eps)
         call f(t,y,ydot) !get derivatives
         do k=1,neq
            yscal(k)=abs(y(k))+abs(h0*ydot(k))+1.d-30
         enddo
         tsec=0.0d0
         call bsstep(y,ydot,neq,tsec,h0,tol,yscal,hdid,hnext,f)
         t=t+tsec/day*TU
         if(abs(h0-hdid).gt.2*eps)then
            !write(0,502) t,h0,hdid,hnext
            h0=day*(tout-t)/TU
            !write(0,*) "new h0",h0
         endif
         !write(0,502) t,h0,hdid,hnext
         !read(5,*)
      enddo

      !write(0,*) j,t,tout
      !write(0,*) h0,hdid,hnext

      tcalc(j)=t  !calculated times for integration model.  This is used in octiming
      do k=1,nbodies
         xpos(k,j)=y(6*k-5) !X
         ypos(k,j)=y(6*k-4) !Y
         zpos(k,j)=y(6*k-3) !Z
      enddo

   enddo
   !read(5,*)


   !for generating plots
   if(iplot.eq.1)then
      do k=1,nbodies
         rx=real(xpos(k,5)-xpos(1,5))
         ry=real(ypos(k,5)-ypos(1,5))
         call pgpt1(rx,ry,-1)
      enddo
   endif

   !if you want the model TTVs, this routine will calculate them.
   call octiming(nbodies,nintg,xpos,ypos,zpos,sol,opos,tc,tcalc,told,bold)


   call transitmodel(nbodies,nintg,xpos,ypos,zpos,sol,tmodel)
   ans(i)=tmodel

   !write(0,*) time(i),ans(i)

   !do k=1,nbodies
   !   write(0,502) m(k)*MU/Mearth,(y(j)*LU/AU,j=6*k-5,6*k-3),(y(j)*LU/TU,j=6*k-5+3,6*k)
      502  format(28(1PE15.8,1X))
   !enddo
   !write(0,*) "Integration loop.."
   !read(5,*)

enddo

if(iplot.eq.1) call pgclos()

return
end
