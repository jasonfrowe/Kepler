subroutine calcbmin()
use precision
implicit none
!import vars

!local vars

!lmdif vars
integer :: lwa,info,n,npt1
integer, allocatable, dimension(:) :: iwa
real(double) :: tollm
real(double), allocatable, dimension(:) :: fvec,wa,bsol
external fcnb


n=1 !number of variables we are fitting 
allocate(bsol(n))
bsol(1)=0.0 !this is the time step

npt1=1 !we only care about impact parameter
lwa=npt1*n+5*npt1*n
allocate(fvec(npt1),iwa(n),wa(lwa))
tollm=1.0d-8
write(0,*) "Calling lmdif1"  !we want to minimize 'b' which is contained in fvec
call lmdif1(fcnb,npt1,n,bsol,fvec,tollm,info,iwa,wa,lwa)
write(0,*) "info: ",info,bsol(1),fvec(1)

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
real(double) :: h0,hdid,eps
!local mercury vars
integer :: dtflag,nclo
integer, parameter :: cmax=50
integer, dimension(cmax) :: iclo,jclo
real(double), dimension(cmax) :: dclo,tclo
real(double), dimension(6,cmax) :: ixvclo,jxvclo


!definition of small number
eps=EPSILON(1.0d0)

!x gives time-step.
h0=x
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
!        if(iplot.eq.1) call pgclos()
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

do j=1,nbodies
    xpos(j,1)=y(6*j-5) !X
    ypos(j,1)=y(6*j-4) !Y
    zpos(j,1)=y(6*j-3) !Z
enddo

call calcimpact(nbodies,y,sol,b_cur)



fvec(1)=bb

return
end