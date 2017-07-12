program koi314rho
use precision
implicit none
integer :: nplot,i,nplotd2,j,nrho,k,nprho,kk
integer, allocatable, dimension(:) :: nprhop
real :: drj,radius,rhostar
real, allocatable, dimension(:) :: rj,x,y
real, allocatable, dimension(:,:) :: rhomodel

nplot=1000
nrho=5
allocate(rhomodel(nrho,3),rj(4),x(nplot),y(nplot))

!Uniform
rhomodel(1,1)=9.5
rhomodel(1,2)=2.0
rhomodel(1,3)=-2.4
!Rayleigh e=0.1
rhomodel(2,1)=9.2
rhomodel(2,2)=1.8
rhomodel(2,3)=-2.2
!Rayleigh e=0.02
rhomodel(3,1)=8.7
rhomodel(3,2)=1.3
rhomodel(3,3)=-1.2
!e=0
rhomodel(4,1)=8.6
rhomodel(4,2)=1.2
rhomodel(4,3)=-1.1
!no priors
rhomodel(5,1)=7.4
rhomodel(5,2)=2.3
rhomodel(5,3)=-3.8

call pgopen('?')
call pgpage()
call PGPAP ( 8.0 ,1.0)
call pgsch(1.2)
call pgslw(3)
call pgvport(0.15,0.85,0.15,0.85)

rj(1)=0.3
rj(2)=0.7
rj(3)=0.3
rj(4)=0.7

call pgwindow(rj(1),rj(2),rj(3),rj(4))
!call pgbox("BCNTS1",0.0,0,"BCNTS1",0.0,0)
call pglabel("Mass (M\d\(2281)\u)","Radius (R\d\(2281)\u)","")
call pgsch(0.8)

call PGSFS(1) !filling style

nprho=2
allocate(nprhop(nprho))
nprhop(1)=1
nprhop(2)=4
do kk=1,nprho
   k=nprhop(kk)
   nplotd2=nplot/2
   drj=(rj(2)-rj(1))/real(nplotd2-1)
   rhostar=rhomodel(k,1)+rhomodel(k,2)
   do i=1,nplotd2
      x(i)=rj(1)+real(i-1)*drj
      y(i)=radius(x(i),rhostar)
   enddo
   rhostar=rhomodel(k,1)+rhomodel(k,3)
   j=-1
   do i=nplotd2+1,nplot
      j=j+1
      x(i)=x(nplotd2-j)
      y(i)=radius(x(i),rhostar)
   enddo
   call pgsci(4+kk)
   call pgpoly(nplot,x,y)
   call pgsci(1)
   call pgline(nplot,x,y)

!   drj=(rj(2)-rj(1))/real(nplot-1)
!   do i=1,nplot
!      x(i)=rj(1)+real(i-1)*drj
!      y(i)=radius(x(i),rhomodel(k,1))
!   enddo
!   call pgline(nplot,x,y)
enddo
call pgsci(1)

call pgsci(15)
!Muirhead  0.51 0.48 0.06 0.06
call pgerr1(5,0.51,0.48,0.06,1.0)
call pgerr1(6,0.51,0.48,0.06,1.0)
!Pineda    0.57 0.54 0.05 0.05
call pgerr1(5,0.57,0.54,0.05,1.0)
call pgerr1(6,0.57,0.54,0.05,1.0)
call pgsci(1)
call pgerr1(5,0.521,0.442,0.055,1.0)
call pgerr1(6,0.521,0.442,0.024,1.0)

call pgsch(1.2)
call pgbox("BCNTS1",0.0,0,"BCNTS1",0.0,0)

call pgclos()

end program koi314rho

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
real function radius(mass,density)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
real :: mass,density,Pi

Pi=acos(-1.d0)!define Pi

radius=(Msun*mass/(4.0/3.0*Pi*density*1000))**(1.0/3.0)/Rsun

return
end

