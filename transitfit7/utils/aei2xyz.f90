subroutine aei2xyz(nbodies,sol,y,m,epoch,percor)
use precision
implicit none
!import vars
integer :: nbodies
real(double) :: epoch
real(double), dimension(:) :: sol,y,m,percor
!local vars
integer :: i,j,np,ii,jj,maxiter,iter
integer, allocatable, dimension(:) :: iPer
real(double) :: adrs,Per,b,Pid2,esinw,ecosw,Psec,varpi,bige,      &
 Meanlong,Meananom,T0,eccanom,lame,bigm,de,eps,bigen,fourpisq,fourpid3,    &
 Rstar,Mstar,rhostar,sqecosw,sqesinw
real(double), allocatable, dimension(:) :: mtot,f,a,ecc,irad,w,Omrad,r, &
 x0,y0,z0,vx0,vy0,vz0,x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,x3,y3,  &
 z3,vx3,vy3,vz3,Pers,x4,y4,z4,vx4,vy4,vz4

eps=EPSILON(1.0d0)

Pid2=Pi/2.0d0
fourpisq=4.0d0*Pi*Pi
fourpid3=4.0d0*Pi/3.0d0

!need to sort by Period...
allocate(Pers(nbodies),iPer(nbodies))
do i=1,nbodies
   np=7+7*(i-1)
   Pers(i)=sol(np+2)
enddo
call rqsort(nbodies,Pers,iPer)
np=7+7*(iPer(1)-1)
rhostar=abs(sol(1)*1000.0) !mean stellar density (kg/m^3)
Mstar=abs(sol(np+5)) !mass of central object (MEarth)
Rstar=(Mearth*Mstar/(fourpid3*rhostar))**(1.0d0/3.0d0)

allocate(mtot(nbodies),a(nbodies),ecc(nbodies),f(nbodies),irad(nbodies),&
 w(nbodies),Omrad(nbodies))
do ii=1,nbodies
   i=iPer(ii) !we work with smallest periods first.
   np=7+7*(i-1)
   Per=sol(np+2)-percor(i) !orbital period in days
   Psec=Per*86400.0d0 !orbital period in sec
   b=sol(np+3) !impact parameter (does not need to be positive)
   sqecosw=sol(np+6)  !e^1/2 cos(w)
   sqesinw=sol(np+7)  !e^1/2 sin(w)
   !make sure eccentricity results in a bounded orbit
   sqecosw=min(0.999,max(-0.999,sqecosw))
   sqesinw=min(0.999,max(-0.999,sqesinw)) 
   ecc(i)=(sqecosw*sqecosw+sqesinw*sqesinw) !eccentricity
   ecosw=sqrt(ecc(i))*sqecosw
   esinw=sqrt(ecc(i))*sqesinw
   a(i)=0.0d0 !initalization of semi-major axis

   Omrad(i)=0.0d0 !should this matter?

   if(ecosw.ne.0.0d0) w(i) = atan(esinw/ecosw)
   if(ecosw.eq.0.0d0) w(i) = asin(esinw)
   if(ecosw.lt.0.0d0) w(i) = w(i)+Pi

!  mass (kg)
   m(i)=abs(sol(np+5)*Mearth)
   if(i.eq.1)then  !star is always first as P=0
      mtot(i)=m(i)
   else
      mtot(i)=m(i)+mtot(iPer(ii-1))
   endif
!   write(0,*) i,mtot(i)

   !calculate inclination
   if(Per>0.0d0)then
      adrs=((Gr*(mtot(i))*Psec*Psec/fourpisq)**(1.0d0/3.0d0))/Rstar
      irad(i)=pid2-acos(b/adrs)
   else
      adrs=0.0d0    !central star does not move
      irad(i)=0.0d0!Pid2
   endif


   if(i.eq.1)then
      a(i)=0.0d0 !star is at center and does not move
      f(i)=0.0d0
   else
      a(i) = (Gr*mtot(i)*Psec*Psec/(4.0*Pi*Pi))**(1.0/3.0) !semi-major axis
      varpi = Omrad(i)+w(i)

      bige = 2.0*atan(sqrt((1.0-ecc(i))/(1.0+ecc(i)))*tan(0.5*(0.5*pi-w(i)))) ! Murray and Dermott eqn 2.46
      Meanlong = varpi+bige-ecc(i)*sin(bige) ! eqn 2.52
      Meananom = Meanlong - varpi         ! M&D 2.53

      T0=sol(np+1) !time of first transit
      eccanom = Meananom+2.0*pi*(epoch-T0)/Per
      lame = eccanom+varpi
      do while (lame.lt.-1.0*pi)
         lame = lame+2.0*pi
      enddo
      do while (lame.gt. 1.0*pi)
         lame = lame-2.0*pi
      enddo

      bigm = lame-varpi
      bige = bigm

      !write(0,*) 'kepler'
      de=1.0d30
      maxiter=250
      iter=0
      do while((de.gt.eps).and.(iter.lt.maxiter)) !iterate until we reach numerical precision
         !write(0,*) 'iter',iter
         bigen = bigm+ecc(i)*sin(bige) !Murray and Dermott eqn 2.54, iterative solution
         de=abs(bigen-bige)
         bige=bigen
         iter=iter+1
      enddo
      !write(0,*) 'kepler done..'
!     do j=1,25
!         bige = bigm+ecc(i)*sin(bige) !Murray and Dermott eqn 2.54, iterative solution
!      enddo

      f(i)=2.0*atan(sqrt((1.0+ecc(i))/(1.0-ecc(i)))*tan(0.5*bige))

   endif

enddo

allocate(r(nbodies),x0(nbodies),y0(nbodies),z0(nbodies),vx0(nbodies),   &
 vy0(nbodies),vz0(nbodies),x1(nbodies),y1(nbodies),z1(nbodies),         &
 vx1(nbodies),vy1(nbodies),vz1(nbodies),x2(nbodies),y2(nbodies),        &
 z2(nbodies),vx2(nbodies),vy2(nbodies),vz2(nbodies),x3(nbodies),        &
 y3(nbodies),z3(nbodies),vx3(nbodies),vy3(nbodies),vz3(nbodies))

y=0 !initialize all values to zero.

do ii=1,nbodies
   i=iPer(ii)
   if(i.eq.1)then
      r(i)=0.0d0  !set all values for central star to zero.
      x0(i)=0.0d0
      y0(i)=0.0d0
      z0(i)=0.0d0
!
      vx0(i)=0.0d0
      vy0(i)=0.0d0
      vz0(i)=0.0d0
!
      x1(i)=0.0d0
      y1(i)=0.0d0
      z1(i)=0.0d0
!
      vx1(i)=0.0d0
      vy1(i)=0.0d0
      vz1(i)=0.0d0
!
      x2(i)=0.0d0
      y2(i)=0.0d0
      z2(i)=0.0d0
!
      vx2(i)=0.0d0
      vy2(i)=0.0d0
      vz2(i)=0.0d0
!
      x3(i)=0.0d0
      y3(i)=0.0d0
      z3(i)=0.0d0
!
      vx3(i)=0.0d0
      vy3(i)=0.0d0
      vz3(i)=0.0d0
!
!      do j=6*i-5,6*i
!         y(j)=0.0d0  !contains x,y,z,vx,vy,vz
!      enddo
   else
      !begin in the orbit plane, 2-d only
      r(i)=a(i)*(1.0-ecc(i)*ecc(i))/(1.0+ecc(i)*cos(f(i))) !M+D eqn 2.20
      x0(i) = r(i)*cos(f(i))
      y0(i) = r(i)*sin(f(i))
      z0(i) = 0.0

      vx0(i) = -sqrt(Gr*(mtot(i))/(a(i)*(1.0-ecc(i)*ecc(i))))*sin(f(i))! 2.36
      vy0(i) = sqrt(Gr*(mtot(i))/(a(i)*(1.0-ecc(i)*ecc(i))))*(ecc(i)+cos(f(i)))
      vz0(i) = 0.0

      !rotate in plane due to omega;
      x1(i) = x0(i)*cos(w(i))- sin(w(i))*y0(i) !P1, eqn 2.119
      y1(i) = x0(i)*sin(w(i))+cos(w(i))*y0(i)
      z1(i) = z0(i)

      vx1(i) = vx0(i)*cos(w(i))- sin(w(i))*vy0(i)
      vy1(i) = vx0(i)*sin(w(i))+cos(w(i))*vy0(i)
      vz1(i) = vz0(i)

      !rotate out of sky-plane due to inclination; P2, eqn 2.119
      x2(i) = x1(i) !the long axis of the ellipse stays the same
      y2(i) = y1(i)*cos(irad(i))-z1(i)*sin(irad(i)) !the y extent of the ellipse shrinks with rotation.
      z2(i) = y1(i)*sin(irad(i))+z1(i)*cos(irad(i)) !the z extent should grow with rotation

      vx2(i) = vx1(i)
      vy2(i) = vy1(i)*cos(irad(i))-vz1(i)*sin(irad(i))
      vz2(i) = vy1(i)*sin(irad(i))+vz1(i)*cos(irad(i))

      !'Precess' plane around reference axes, P3, eqn 2.121, make negative as in line 74 setup.pro
      x3(i) = -1.0*((x2(i)*cos(Omrad(i)) - y2(i)*sin(Omrad(i))))
      y3(i) = -1.0*(x2(i)*sin(Omrad(i)) + y2(i)*cos(Omrad(i)))
      z3(i) = -1.0*(z2(i))

      vx3(i) = -1.0*(vx2(i)*cos(Omrad(i)) - vy2(i)*sin(Omrad(i)))
      vy3(i) = -1.0*(vx2(i)*sin(Omrad(i)) + vy2(i)*cos(Omrad(i)))
      vz3(i) = -1.0*(vz2(i))

   endif

enddo

allocate(x4(nbodies),y4(nbodies),z4(nbodies),vx4(nbodies),vy4(nbodies), &
 vz4(nbodies))

do i=1,nbodies
   x4(i)=0.0
   y4(i)=0.0
   z4(i)=0.0
   vx4(i)=0.0
   vy4(i)=0.0
   vz4(i)=0.0
enddo

do ii=1,nbodies
   i=iPer(ii)
   x4(i)=x3(i)
   y4(i)=y3(i)
   z4(i)=z3(i)
   vx4(i)=vx3(i)
   vy4(i)=vy3(i)
   vz4(i)=vz3(i)
   do jj=1,ii-1
      j=iPer(jj)
      x4(i)=x4(i)+x3(j)*m(j)/mtot(j)
      y4(i)=y4(i)+y3(j)*m(j)/mtot(j)
      z4(i)=z4(i)+z3(j)*m(j)/mtot(j)
      vx4(i)=vx4(i)+vx3(j)*m(j)/mtot(j)
      vy4(i)=vy4(i)+vy3(j)*m(j)/mtot(j)
      vz4(i)=vz4(i)+vz3(j)*m(j)/mtot(j)
   enddo
   y(6*i-5)=x4(i)
   y(6*i-4)=-y4(i)
   y(6*i-3)=-z4(i)
   y(6*i-2)=vx4(i)
   y(6*i-1)=-vy4(i)
   y(6*i  )=-vz4(i)
enddo

!write(0,*) "nbodies: ",nbodies
do i=1,nbodies
!   write(0,*) "MU: ",MU
   m(i)=m(i)/MU*K2
   y(6*i-5)=y(6*i-5)/LU
   y(6*i-4)=y(6*i-4)/LU
   y(6*i-3)=y(6*i-3)/LU
   y(6*i-2)=y(6*i-2)*TU/LU
   y(6*i-1)=y(6*i-1)*TU/LU
   y(6*i)=y(6*i)*TU/LU
   !write(0,500) m(i)/K2*Msun/Mearth,(y(j),j=6*i-5,6*i),a(i)/LU
   500  format(28(1PE10.3,1X))
enddo

end subroutine aei2xyz
