subroutine aei2xyz(nbodies,sol,y,epoch)
use precision
implicit none
!import vars
integer :: nbodies
real(double) :: epoch
real(double), dimension(:) :: sol,y
!local vars
integer :: i,j,k,np,ii,jj
integer, allocatable, dimension(:) :: iPer
real(double) :: eps,Pi,Pid2,tPi,fourpisq,fourpid3,Per,Psec,b,ecosw,esinw,t0,Eanom,phi0,eccn,w,t,phi,Manom,Tanom, &
 trueanomaly,Mtotal,Mp,u,at,Ms,ap,as,rhostar,Mstar,Rstar,distance,rp,rs,astar,drs,incl,dumr,Gr,vel
real(double), allocatable, dimension(:) :: Pers
real(double), allocatable, dimension(:,:) :: pcart,scart,pcart2,scart2

eps=EPSILON(1.0d0)

Pi=3.141592653589793d0!define Pi
tPi=2.0d0*Pi
Pid2=Pi/2.0d0
fourpisq=4.0d0*Pi*Pi
fourpid3=4.0d0*Pi/3.0d0
Gr=6.6719842229637d-11

!need to sort by Period...
allocate(Pers(nbodies),iPer(nbodies))
do i=1,nbodies
   np=7+7*(i-1)
   Pers(i)=sol(np+2)
   m(i)=abs(sol(np+5))*Mearth !mass of each object (kg)
enddo
call rqsort(nbodies,Pers,iPer)

np=7+7*(iPer(1)-1)
rhostar=abs(sol(1)*1000.0) !mean stellar density (kg/m^3)
Mstar=abs(sol(np+5)) !mass of central object (MEarth)
Rstar=(Mearth*Mstar/(fourpid3*rhostar))**(1.0d0/3.0d0)
Mtotal=abs(sol(np+5)) !initiate sum of masses [Mearth]

allocate(pcart(6,nbodies),scart(6,nbodies))
pcart=0.0d0
scart=0.0d0

do ii=2,nbodies !The central body will have a period of zero. 

   i=iPer(ii) !we work with smallest periods first.

   np=7+7*(i-1)
   Per=sol(np+2) !orbital period in days
   Psec=Per*86400.0d0 !orbital period in sec
   b=sol(np+3) !impact parameter
   ecosw=sol(np+6)  !e cos(w)
   esinw=sol(np+7)  !e sin(w)

   eccn=sqrt(ecosw*ecosw+esinw*esinw) !eccentricity
   if(ecosw.ne.0.0d0) w = atan(esinw/ecosw)
   if(ecosw.eq.0.0d0) w = asin(esinw)
   if(ecosw.lt.0.0d0) w = w+Pi

   t0=sol(np+1) !mid-time of transit
   Eanom=tan(w/2.0d0)/sqrt((1.0d0+eccn)/(1.0d0-eccn)) !mean anomaly
   Eanom=2.0d0*atan(Eanom)
   phi0=Eanom-eccn*sin(Eanom) 

   t=epoch-t0 !current time (epoch) relative to centre of transit time (t0)

   !get orbital position (mean anomaly)
   phi=t/per-floor(t/per)
   phi=phi*tPi+phi0
   Manom=phi
   if(Manom.gt.tPi) Manom=Manom-tPi
   if(Manom.lt.0.0d0) Manom=Manom+tPi
   !get True anomaly 
   call kepler(Manom,Eanom,eccn)
   Tanom=trueanomaly(eccn,Eanom)

   !get semi-major axis from Kepler's Third Law.
   Mp=abs(sol(np+5))  !mass of planet (Mearth)
   Ms=Mtotal !enclosed mass of central object (Mearth)
   Mtotal=Mtotal+Mp !total mass with new object (Mearth)
   u=MEarth*Gr*Mtotal !G(m1+m2) [MKS] 
   at=u*Psec*Psec/fourpisq 
   at=at**(1.0d0/3.0d0)  !semi-major axis (meters)

   ap = at*Ms/Mtotal !semi-major axis of planet (meters)
   as = at*Mp/Mtotal !semi-mjaor axis of central source (meters)

   rp = distance(ap,eccn,Tanom) !planet distance from CoM
   rs = rp*Mp/Ms !distance(as,eccn,Tanom) !star distance from CoM

   !get inclination of orbit relative to Sun.
   astar=(Mearth*Gr*(Mstar+Mp)*Psec*Psec/fourpisq)**(1.0d0/3.0d0) ![m]
   drs=astar/Rstar
   incl=acos(b/drs) !drs is rp/R*, which we get from rhostar.  Rstar=(Mstar/(4pi*rhostar))**(1/3)
   !write(0,*) "incl ",b,incl

   !so far CoM is located at x,y=0.
   !planet position in CoM 
   pcart(1,i)=rp*sin(Tanom-w)
   pcart(2,i)=rp*Cos(Tanom-w)*sin(incl)
   pcart(3,i)=rp*Cos(Tanom-w)*cos(incl)

   !planet velocities in CoM
   !dumr=tPi*rp/(Psec*sqrt(1-eccn*eccn))
   !pcart(4,i)=dumr*(eccn+cos(Tanom-w))
   !pcart(5,i)=-dumr*sin(Tanom-w)*sin(incl)
   !pcart(6,i)=-dumr*sin(Tanom-w)*cos(incl)
   vel=sqrt(u*(2/rp-1/(rp+rs)))
   pcart(4,i)=vel*cos(Tanom-w)
   pcart(5,i)=-vel*sin(Tanom-w)*sin(incl)
   pcart(6,i)=-vel*sin(Tanom-w)*cos(incl)


   !Star positions in CoM
   do j=1,6
      scart(j,i)=-Mp*pcart(j,i)/Ms
   enddo
   !star velocities in CoM
   !dumr=tPi*rs/(Psec*sqrt(1-eccn*eccn))
   !scart(4,i)=-dumr*(eccn+cos(Tanom-w))
   !scart(5,i)=dumr*sin(Tanom-w)*sin(incl)
   !scart(6,i)=dumr*sin(Tanom-w)*cos(incl)

   write(0,500) i,m(i),Per,ap/AU,as/AU,astar/AU
   write(0,501) (pcart(j,i)/AU,j=1,3),(pcart(j,i),j=4,6)
   write(0,501) (scart(j,i)/AU,j=1,3),(scart(j,i),j=4,6)
   500 format(i1,1x,50(1X,1PE17.10))
   501 format(2x,50(1X,1PE17.10))

enddo


allocate(pcart2(6,nbodies),scart2(6,nbodies))

pcart2=pcart
scart2=scart

do ii=2,nbodies
   i=iPer(ii) !we work with smallest periods first.

   do jj=ii+1,nbodies
      j=iPer(jj)
      !write(0,*) i,j
      do k=1,6
         pcart2(k,i)=pcart2(k,i)+scart(k,j)
         scart2(k,i)=scart2(k,i)+scart(k,j)
      enddo
   enddo

   if(ii.eq.2)then
      y(6*1-5)=scart2(1,i)
      y(6*1-4)=scart2(2,i)
      y(6*1-3)=scart2(3,i)
      y(6*1-2)=scart2(4,i)
      y(6*1-1)=scart2(5,i)
      y(6*1)=scart2(6,i)
   endif

   y(6*i-5)=pcart2(1,i)
   y(6*i-4)=pcart2(2,i)
   y(6*i-3)=pcart2(3,i)
   y(6*i-2)=pcart2(4,i)
   y(6*i-1)=pcart2(5,i)
   y(6*i)=pcart2(6,i)

   !write(0,500) i,m(i)/Mearth,Pers(i),ap/AU,as/AU,astar/AU
   !write(0,501) (pcart2(j,i)/AU,j=1,3),(pcart2(j,i),j=4,6)
   !write(0,501) (scart2(j,i)/AU,j=1,3),(scart2(j,i),j=4,6)
   !500 format(i1,1x,50(1X,1PE17.10))
   !501 format(2x,50(1X,1PE17.10))

enddo


do i=1,nbodies
   m(i)=m(i)/MU*K2
   y(6*i-5)=y(6*i-5)/LU
   y(6*i-4)=y(6*i-4)/LU
   y(6*i-3)=y(6*i-3)/LU
   y(6*i-2)=y(6*i-2)*TU/LU
   y(6*i-1)=y(6*i-1)*TU/LU
   y(6*i)=y(6*i)*TU/LU
   write(0,502) m(i)/K2*Msun/Mearth,(y(j),j=6*i-5,6*i)!,a(i)/LU
   502  format(28(1PE10.3,1X))
enddo

!write(0,*) "Reached End of aei2xyz"
!read(5,*)

return
end subroutine aei2xyz
