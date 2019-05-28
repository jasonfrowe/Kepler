subroutine lcmodel(nbodies,npt,tol,sol,time,itime,ntmid,tmid,percor,ans,colflag,itprint,itmodel)
use precision
implicit none
!import vars
integer :: nbodies,npt,itprint,itmodel,colflag
real(double) :: tol
real(double), dimension(:) :: sol,time,itime,percor,ans
integer, dimension(:) :: ntmid !used with octiming 
real(double), dimension(:,:) :: tmid !used with octiming
!local vars
integer :: nintg,i,j,k,neq,lrw,liw,istate,nbod
integer, allocatable, dimension(:) :: iwork,tc
real(double) :: dnintg,tdnintg,t,tout,tout2,jm1,tmodel,   &
 epoch
real(double), allocatable, dimension(:) :: y,m,rwork,tcalc,opos
real(double), allocatable, dimension(:,:) :: xpos,ypos,zpos,told,bold
!mercury vars
integer :: nclo
integer, parameter :: cmax=50
integer, dimension(cmax) :: iclo,jclo
real(double), dimension(cmax) :: dclo,tclo
real(double), dimension(6,cmax) :: ixvclo,jxvclo
integer :: algor,nbig,ngflag,opflag,dtflag
integer, dimension(8) :: opt
integer, allocatable, dimension(:) :: stat
real(double) :: h0,rcen,rmax,cefac,tstart,hdid,eps,maxint
real(double), dimension(3) :: jcen,en,am
real(double), allocatable, dimension(:) :: rho,rceh,rphys,rce,rcrit
real(double), allocatable, dimension(:,:) :: xh,vh,s,x,v
logical :: first
!plotting stuff
integer :: iplot
real :: bb(4),rasemi,rx,ry
external f,jac
!timing output
integer :: nunit,ip,np
real(double) :: ttomc,t0,per
character(80) :: filename
!save vars for Mercury
integer, parameter :: nmax=2000
real(double) :: a(3,nmax),hrec,angf(3,nmax),ausr(3,nmax)


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
    told,bold,ntmid,tmid)
      use precision
      implicit none
      integer, intent(in) :: nbodies,nintg
      integer, dimension(:), intent(inout) :: tc
      integer, dimension(:), intent(inout) :: ntmid
      real(double), dimension(:,:), intent(in) :: xpos,ypos,zpos
      real(double), dimension(:), intent(in) :: sol,tcalc
      real(double), dimension(:), intent(inout) :: opos
      real(double), dimension(:,:), intent(inout) :: told,bold
      real(double), dimension(:,:), intent(inout) :: tmid
   end subroutine octiming
end interface

!plotting
iplot=0 !set iplot=1 for plotting 

!for outputting timing measurements
!itprint=1  !0 - do not output to file, 1 - output to 'tmeas.dat'

!definition of small number
eps=EPSILON(1.0d0)

!inititation for timing measurements
ntmid=0

!oversampling
nintg=11
dnintg=dble(nintg)
tdnintg=2.0d0*dnintg

allocate(opos(nbodies),tc(nbodies)) !used by octiming for previous position
allocate(tcalc(nintg)) !storing model times for octiming
allocate(told(nbodies,5),bold(nbodies,5))
opos=0.0d0 !initilization to zero
tc=0 !initializtion of transit conidition flag to zero.

!allocate array to contain masses and position,velocity for each body.
allocate(y(nbodies*6),m(nbodies)) !need an extra space to avoid segs (?!)
!do i=1,nbodies*6
!   y(i)=0.0d0
!enddo
y=0.0d0
m=0.0d0
epoch=time(1)
!convert from Keplerian to Cartesian coords 
call aei2xyz(nbodies,sol,y,m,epoch,percor)
!read(5,*)

!!plotting stuff...
!if(iplot.eq.1)then 
!   call pgopen('/xserve') !open PGPlot device
!   call pgpap(8.0,1.0) !page size
!   call pgslw(3) !thicker lines
!   rasemi=0.0
!   do k=1,nbodies
!      rasemi=max(rasemi,real(abs(y(6*k-5)))) !X
!      rasemi=max(rasemi,real(abs(y(6*k-4)))) !Y
!      rasemi=max(rasemi,real(abs(y(6*k-3)))) !Z
!   enddo
!   bb(1)=-3.0*rasemi
!   bb(2)= 3.0*rasemi
!   bb(3)=-3.0*rasemi
!   bb(4)= 3.0*rasemi
!   call pgwindow(bb(1),bb(2),bb(3),bb(4))
!   CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
!   call pgslw(1) !thicker lines
!endif

!arrays to contaim time stamps and x,y,z positions of the bodies
allocate(xpos(nbodies,nintg),ypos(nbodies,nintg),zpos(nbodies,nintg))

nbod=nbodies
neq=6*nbodies !number of equations 6 times number of particles
t=time(1)!/365.25!/TU*86400.0d0     !time zero

!Mercury vars
algor=10 !algorithm being used.
jcen=0.0d0 !J2,J4,J6 of central body
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
rceh=1.0d0 !close-encounter limit (Hill radii)
ngflag=0 !non-grav forces is turned off.
opt=0 !options
en=0 !energies
am=0 !angular momentums
stat=0 !0-alive,1-to be removed
opflag=0 !integration mode
colflag=0 !collision flag.  0 means no collision

maxint=maxintg/86400.0 !sampling [days]  !1-5 min seems to be fine for Kepler.


ans=0.0 !default, set answer array to 0 flux.
!write(0,*) "Starting mercury.."
first=.true.
do i=1,npt

   !target integration time.. (first step is negative!)
   tout=(time(i)-itime(i)*(0.5d0-1.0d0/tdnintg))!/365.25!/TU*86400.0d0

!   write(0,*) 'tout:',t,tout
!   write(6,501) t,(y(6*k-5),y(6*k-4),y(6*k-3),k=1,nbodies)

   h0=tout-t

   if(first) then !setting up co-ordinates and close encounters
      first=.false.
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
      tstart=t !the value of tstart does not matter.. 
   else
      dtflag=2 !normal call.
   endif


   hdid=0.0d0
   do while(abs(h0).gt.2.0*eps)
      hdid=min(maxint,h0)
      call mdt_hy (t,tstart,hdid,tol,rmax,en,am,jcen,rcen,nbod,           &
         nbig,m,x,v,s,rphys,rcrit,rce,stat,algor,opt,dtflag,        &
         ngflag,opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo, &
         a,hrec,angf,ausr)
      if(colflag.ne.0)then
         write(0,*) "Close encounter, stopping integration"
!         if(iplot.eq.1) call pgclos()
         return
      endif
      t=t+hdid
      h0=h0-hdid
   enddo

   call mco_dh2h (t,jcen,nbod,nbig,h0,m,x,v,xh,vh)
   do j=1,nbodies
      y(6*j-5)=xh(1,j)
      y(6*j-4)=xh(2,j)
      y(6*j-3)=xh(3,j)
      y(6*j-2)=vh(1,j)
      y(6*j-1)=vh(2,j)
      y(6*j)=vh(3,j)
   enddo

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

      hdid=0.0d0
      do while(abs(h0).gt.2.0*eps)
         hdid=min(maxint,h0)
         dtflag=2 !normal call.
         call mdt_hy (t,tstart,hdid,tol,rmax,en,am,jcen,rcen,nbod,           &
            nbig,m,x,v,s,rphys,rcrit,rce,stat,algor,opt,dtflag,        &
            ngflag,opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo, &
            a,hrec,angf,ausr)
         if(colflag.ne.0)then
            write(0,*) "Close encounter, stopping integration"
!            if(iplot.eq.1) call pgclos()
            return
         endif
         t=t+hdid
         h0=h0-hdid
      enddo

      call mco_dh2h (t,jcen,nbod,nbig,h0,m,x,v,xh,vh) !t value is not needed 
      do k=1,nbodies
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
!   write(6,501) time(i),(xpos(k,2),ypos(k,2),zpos(k,2),k=1,nbodies)
   501  format(78(1PE13.6,1X))
!   read(5,*)

!   !for generating plots
!   if(iplot.eq.1)then
!      do k=1,nbodies
!         rx=real(xpos(k,5)-xpos(1,5))
!         ry=real(ypos(k,5)-ypos(1,5))
!         call pgpt1(rx,ry,-1)
!      enddo
!   endif

   !if you want the model TTVs, this routine will calculate them.
   call octiming(nbodies,nintg,xpos,ypos,zpos,sol,opos,tc,tcalc,told,bold,ntmid,tmid)

   if(itmodel.eq.1)then
      call transitmodel(nbodies,nintg,xpos,ypos,zpos,sol,tmodel)
      ans(i)=tmodel
   endif

enddo
!write(0,*) "Finished mercury"

!update this to output timing for all nbodies with per>0
if(itprint.eq.1)then
   nunit=12
   do ip=1,nbodies
      np=7+7*(ip-1)
      t0=sol(np+1)
      Per=sol(np+2)
      if((per.gt.0.0d0).and.(ntmid(ip).gt.0))then
         write(filename,'(A6,I1,A4)') 'tmeas_',ip,'.dat'
         open(unit=nunit,file=filename) !open file if we are outputting tt measurements
         do i=1,ntmid(ip)
            ttomc=(tmid(ip,i)-t0)/per-int((tmid(ip,i)-t0)/per)
            if (ttomc.gt.0.5) ttomc=ttomc-1.0
            ttomc=Per*ttomc !convert from phase to days
            write(nunit,502) tmid(ip,i),ttomc
         enddo
   502   format(78(1PE13.6,1X))
         close(nunit)
      endif
   enddo
endif

!if(iplot.eq.1) call pgclos()

return
end
