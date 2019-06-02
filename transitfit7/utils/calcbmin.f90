subroutine calcbmin(nbc,nbodies,t,sol,tol,nbod,m,x,v,algor,nbig,ngflag,opflag,colflag,opt,stat,rcen,rmax,&
    tstart,jcen,en,am,rphys,rce,rcrit,s,a,hrec,angf,ausr,ttran)
use precision
use bfittingmod
implicit none
!import vars
integer, target :: nbc,nbodies
real(double), target :: t
real(double) :: ttran
real(double), dimension(:), target :: sol
!import mercury vars
integer, target :: algor,nbig,ngflag,opflag,colflag,nbod
integer, dimension(8), target :: opt
integer, dimension(:), target :: stat
real(double), target :: rcen,rmax,tstart,tol
real(double), dimension(3), target :: jcen,en,am
real(double), dimension(:), target :: rphys,rce,rcrit,m
real(double), dimension(:,:), target :: s,x,v
!import mercury save vars
integer, parameter :: nmermax=2000
real(double), target :: a(3,nmermax),hrec,angf(3,nmermax),ausr(3,nmermax)
!local vars

!lmdif vars
integer :: lwa,info,n,npt1
integer, allocatable, dimension(:) :: iwa
real(double) :: tollm
real(double), allocatable, dimension(:) :: fvec,wa,bsol
external fcn


n=1 !number of variables we are fitting 
allocate(bsol(n))
bsol(1)=0.0 !this is the time step

!pointers for all of the Mercury vars
nbc2=>nbc
nbodies2=>nbodies
t2=>t
sol2=>sol
tol2=>tol
algor2=>algor
nbig2=>nbig
ngflag2=>ngflag
opflag2=>opflag
colflag2=>colflag
nbod2=>nbod
opt2=>opt
stat2=>stat
rcen2=>rcen
rmax2=>rmax
tstart2=>tstart
jcen2=>jcen
en2=>en
am2=>am
rphys2=>rphys
rce2=>rce
rcrit2=>rcrit
m2=>m
s2=>s
x2=>x
v2=>v
a2=>a
hrec2=>hrec
angf2=>angf
ausr2=>ausr

npt1=1 !we only care about impact parameter
lwa=npt1*n+5*npt1+n
allocate(fvec(npt1),iwa(n),wa(lwa))
tollm=1.0d-8
!write(0,*) "Calling lmdif1"  !we want to minimize 'b' which is contained in fvec
call lmdif1(fcn,npt1,n,bsol,fvec,tollm,info,iwa,wa,lwa)
!write(0,*) "info: ",info,bsol(1),fvec(1)

ttran=bsol(1)+t

nullify(nbc2,nbodies2,t2,sol2,tol2,algor2,nbig2,ngflag2,opflag2,colflag2,nbod2,opt2,stat2,&
  rcen2,rmax2,tstart2,jcen2,en2,am2,rphys2,rce2,rcrit2,m2,s2,x2,v2,a2,hrec2,angf2,ausr2)

deallocate(bsol,fvec,iwa,wa)

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn(npt,n,x,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use bfittingmod
implicit none
!import vars
integer :: npt,n,iflag
real(double), dimension(n) :: x
real(double), dimension(npt) :: fvec
!local vars
integer :: j
real(double) :: h0,hdid,eps,maxint
real(double), allocatable, dimension(:) :: b_cur
!local mercury vars
integer :: dtflag,nclo
integer, parameter :: cmax=50
integer, dimension(cmax) :: iclo,jclo
real(double), dimension(cmax) :: dclo,tclo
real(double), dimension(6,cmax) :: ixvclo,jxvclo
!local mercury grabs
real(double) :: tb
real(double), allocatable, dimension(:) :: m3,y
real(double), allocatable, dimension(:,:) :: xh3,vh3,x3,v3
integer, parameter :: nmermax=2000
real(double), target :: a3(3,nmermax),hrec3,angf3(3,nmermax),ausr3(3,nmermax)

interface
   subroutine calcimpact(nbodies,y,sol,b_cur)
      use precision
      integer :: nbodies
      real(double), dimension(:) :: y,sol,b_cur
   end subroutine calcimpact
end interface

!grab pointers
tb=t2
allocate(m3(nbodies2),xh3(3,nbod2),vh3(3,nbod2),x3(3,nbod2),v3(3,nbod2),y(nbodies2*6))
y=0
m3=m2
x3=x2
v3=v2
a3=a2
hrec3=hrec2
angf3=angf2
ausr3=ausr2

allocate(b_cur(nbodies2))

!definition of small number
eps=EPSILON(1.0d0)
maxint=maxintg/86400.0 !default for out of transit sampling [days] 

!x gives time-step.
h0=x(1)
hdid=0.0d0
do while(abs(h0).gt.2.0*eps)
 	hdid=min(maxint,h0)
    dtflag=2 !normal call.
    call mdt_hy (tb,tstart2,hdid,tol2,rmax2,en2,am2,jcen2,rcen2,nbod2,           &
     nbig2,m3,x3,v3,s2,rphys2,rcrit2,rce2,stat2,algor2,opt2,dtflag,        &
     ngflag2,opflag2,colflag2,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo, &
     a3,hrec3,angf3,ausr3)
    if(colflag2.ne.0)then
        write(0,*) "Close encounter, stopping integration"
!        if(iplot.eq.1) call pgclos()
        return
    endif
    tb=tb+hdid
    h0=h0-hdid
enddo
call mco_dh2h (tb,jcen2,nbod2,nbig2,h0,m3,x3,v3,xh3,vh3)
do j=2,nbodies2
    y(6*j-5)=xh3(1,j)
    y(6*j-4)=xh3(2,j)
    y(6*j-3)=xh3(3,j)
    y(6*j-2)=vh3(1,j)
    y(6*j-1)=vh3(2,j)
    y(6*j)  =vh3(3,j)
enddo

!write(0,*) y
call calcimpact(nbodies2,y,sol2,b_cur)

fvec(1)=b_cur(nbc2)

!write(0,*) "bb: ",b_cur
!write(0,*) fvec(1),x(1)
!read(5,*)

deallocate(m3,xh3,vh3,x3,v3,y,b_cur)

return
end
