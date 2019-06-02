subroutine lcmodel_pc(nbodies,npt,tol,sol,time,ntmid,tmid,percor,colflag,itprint)
use precision
implicit none
!import vars
integer :: nbodies,colflag,itprint,npt
real(double) :: tol
real(double), dimension(:) :: sol,time,percor
integer, dimension(:) :: ntmid !used with octiming 
real(double), dimension(:,:) :: tmid !used with octiming
!local vars
integer :: nintg,nbod,neq,i,j
real(double) :: eps,epoch,t,te,tout
integer, allocatable, dimension(:) :: tc
real(double), allocatable, dimension(:) :: tcalc,opos,y,m
real(double), allocatable, dimension(:,:) :: told,bold,xpos,ypos,zpos
!mercury vars
integer :: nclo
integer, parameter :: cmax=50
integer, dimension(cmax) :: iclo,jclo
real(double), dimension(cmax) :: dclo,tclo
real(double), dimension(6,cmax) :: ixvclo,jxvclo
integer :: algor,nbig,ngflag,opflag,dtflag
integer, dimension(8) :: opt
integer, allocatable, dimension(:) :: stat
real(double) :: h0,rcen,rmax,cefac,tstart,maxint
real(double), dimension(3) :: jcen,en,am
real(double), allocatable, dimension(:) :: rho,rceh,rphys,rce,rcrit
real(double), allocatable, dimension(:,:) :: xh,vh,s,x,v
logical :: first
!!plotting stuff
integer :: iplot
!integer :: k
!real :: bb(4),rasemi,rx,ry
!vars to control integration steps
real(double) :: RpRs
real(double), allocatable, dimension(:) :: b_thres,b_old,b_cur
!save vars for Mercury
real(double) :: a(3,nmermax),hrec,angf(3,nmermax),ausr(3,nmermax)
!timing output
integer :: nunit,ip,np
real(double) :: ttomc,t0,per
character(80) :: filename
character(200) :: dumc


interface
   subroutine aei2xyz(nbodies,sol,y,m,epoch,percor)
      use precision
      implicit none
      integer, intent(inout) :: nbodies
      real(double), dimension(:), intent(inout) :: sol,percor
      real(double), dimension(:), intent(inout) :: y,m
      real(double), intent(inout) :: epoch
   end subroutine aei2xyz
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
   subroutine calcimpact(nbodies,y,sol,b_cur)
      use precision
      integer, intent(in) :: nbodies
      real(double), dimension(:), intent(in) :: y,sol
      real(double), dimension(:), intent(inout) ::b_cur
   end subroutine calcimpact
end interface

!plotting
iplot=0 !set iplot=1 for plotting 

!definition of small number
eps=EPSILON(1.0d0)

!inititation for timing measurements
ntmid=0

!oversampling
nintg=1 !no oversampling needed 
!dnintg=dble(nintg)
!tdnintg=2.0d0*dnintg

allocate(opos(nbodies),tc(nbodies)) !used by octiming for previous position
allocate(tcalc(nintg)) !storing model times for octiming
allocate(told(nbodies,5),bold(nbodies,5))
opos=0.0d0 !initilization to zero
tc=0 !initializtion of transit conidition flag to zero.

!allocate array to contain masses and position,velocity for each body.
allocate(y(nbodies*6),m(nbodies))
y=0.0d0
m=0.0d0

epoch=time(1)
!convert from Keplerian to Cartesian coords 
call aei2xyz(nbodies,sol,y,m,epoch,percor)

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
xpos=0.0d0
ypos=0.0d0
zpos=0.0d0

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

maxint=maxintg/86400.0 !default for out of transit sampling [days] 
!check for transit condition and if so, adjust integration time for timing precision
!calculate thresholds for a transit to occur
allocate(b_thres(nbodies),b_old(nbodies),b_cur(nbodies))
do i=2,nbodies
   np=7+7*(i-1)
   RpRs=abs(sol(np+4))
   b_thres(i)=1.0+RpRs
enddo
call calcimpact(nbodies,y,sol,b_cur)

b_old=b_cur !initalize b_old
do i=2,nbodies
   if(b_cur(i).lt.b_thres(i)) maxint=maxintg_nt/86400.0
enddo
!write(0,*) (b_cur(i),i=2,nbodies)
!write(0,*) maxint*86400.0 
!read(5,*)


te=maxval(time(1:npt))

first=.true.
do while(t.le.te)

   tout=t+maxint !target integration time
   h0=tout-t     !target integration step

   if(first) then !setting up co-ordinates and close encounters
       first=.false.
         do j=1,nbodies
            xh(1,j)=y(6*j-5)
            xh(2,j)=y(6*j-4)
            xh(3,j)=y(6*j-3)
            vh(1,j)=y(6*j-2)
            vh(2,j)=y(6*j-1)
            vh(3,j)=y(6*j)
         enddo
         call mce_init (t,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,       &
          m,xh,vh,s,rho,rceh,rphys,rce,rcrit,1)
         call mco_h2dh (t,jcen,nbod,nbig,h0,m,xh,vh,x,v)
         dtflag=0 !first time calling mdt_hy
         tstart=t !the value of tstart does not matter.. 
      else
         dtflag=2 !normal call.
      endif

       call mdt_hy (t,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,           &
        nbig,m,x,v,s,rphys,rcrit,rce,stat,algor,opt,dtflag,        &
        ngflag,opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo, &
        a,hrec,angf,ausr)
    if(colflag.ne.0)then
        write(0,*) "Close encounter, stopping integration"
!        if(iplot.eq.1) call pgclos()
        return
    endif
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

    tcalc(1)=t
    do j=1,nbodies
      xpos(j,1)=y(6*j-5) !X
      ypos(j,1)=y(6*j-4) !Y
      zpos(j,1)=y(6*j-3) !Z
   enddo

!   !for generating plots
!   if(iplot.eq.1)then
!      do k=1,nbodies
!         rx=real(xpos(k,1)-xpos(1,1))
!         ry=real(ypos(k,1)-ypos(1,1))
!         call pgpt1(rx,ry,-1)
!      enddo
!   endif

   call octiming(nbodies,nintg,xpos,ypos,zpos,sol,opos,tc,tcalc,told,bold,&
    ntmid,tmid)


    maxint=maxintg/86400.0 !default for out of transit sampling [days] 
    !check for transit condition and if so, adjust integration time for timing precision
    !calculate thresholds for a transit to occur
    call calcimpact(nbodies,y,sol,b_cur)

    do j=2,nbodies
        !if(b_cur(j).lt.2.0d0)then 
        if(abs(2.0*b_cur(j)-b_old(j)).lt.b_thres(j))then
            maxint=maxintg_nt/86400.0
            !write(0,*) b_cur(j),maxint*86400.0 
            !read(5,*)
        endif
    enddo
    b_old=b_cur !update b_old

    write(dumc,503) 't: ',t,maxint,xpos(2,1),ypos(2,1),zpos(2,1)
 503   format(A3,1X,78(1PE13.6,1X))
!    !read(5,*)

enddo


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
            if(ip.eq.3)then
            	write(6,502) tmid(ip,i),ttomc
            endif
         enddo
   502   format(78(1PE13.6,1X))
         close(nunit)
      endif
   enddo
endif

!if(iplot.eq.1) call pgclos()

return 
end subroutine
