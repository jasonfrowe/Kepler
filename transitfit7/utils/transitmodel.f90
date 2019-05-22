subroutine transitmodel(nbodies,nintg,xpos,ypos,zpos,sol,tmodel)
use precision
implicit none
integer nbodies,nintg,i,j,caltran,np
real(double) :: RpRs,rstar,LU2,tmodel,c1,c2,c3,c4,zpt,fourpid3,Mstar,rhostar
real(double), allocatable, dimension(:) :: b,tflux,mum
real(double), dimension(:,:) :: xpos,ypos,zpos
real(double), dimension(:) :: sol

LU2=LU*LU

allocate(b(nintg),tflux(nintg),mum(nintg))

!rstar=(sol(12)/(4.0d0/3.0d0*Pi*sol(1)*1000.0d0)*Mearth)**(1.0d0/3.0d0)
fourpid3=4.0d0*Pi/3.0d0
rhostar=abs(sol(1)*1000.0) !mean stellar density (kg/m^3)
Mstar=abs(sol(12)) !mass of central object (MEarth)
rstar=(Mearth*Mstar/(fourpid3*rhostar))**(1.0d0/3.0d0)

!write(0,*) "Rstar: ",rstar/Rsun,sol(12),sol(12)*MU/Msun
c1=sol(2)!limb-darkening
c2=sol(3)
c3=sol(4)
c4=sol(5)
zpt=sol(7)!zero point

tmodel=0.0d0
do i=2,nbodies
   np=7+7*(i-1)
!   write(0,*) "test"
   RpRs=abs(sol(np+4))
!   write(0,*) "RpRs: ",RpRs
   if(RpRs.le.0.0d0)then
      tmodel=tmodel+0.0d0
      cycle
   endif

   if(ypos(i,nintg/2)-ypos(1,nintg/2).gt.0.0d0)then !we have a transit
      do j=1,nintg
         b(j)=(zpos(i,j)-zpos(1,j))*(zpos(i,j)-zpos(1,j))+ &
            (xpos(i,j)-xpos(1,j))*(xpos(i,j)-xpos(1,j))
         b(j)=sqrt(LU2*b(j))/rstar
!         write(6,*) j,b(j)
      enddo

      caltran=0
      do j=1,nintg
         if(b(j).lt.1.0d0+RpRs)then
            caltran=1
         endif
      enddo

      if(caltran.eq.1) then
         if((c3.eq.0.0).and.(c4.eq.0.0))then
            call occultquad(b,c1,c2,RpRs,tflux,mum,nintg)
            !Kipping co-efficients
         elseif((c1.eq.0.0).and.(c2.eq.0.0))then
            c1=2.0d0*sqrt(c3)*c4 !convert to regular LD
            c2=sqrt(c3)*(1.0d0-2.0d0*c4)
            call occultquad(b,c1,c2,RpRs,tflux,mum,nintg)
            c1=0.0d0  !zero out entries.
            c2=0.0d0
         else
          !non-linear law.
            call occultsmall(RpRs,c1,c2,c3,c4,nintg,b,tflux)
!            call occultnl(RpRs,c1,c2,c3,c4,b,tflux,
!     .      mulimbf,nintg)
!            write(0,*) RpRs,abs(sol(10*(ii-1)+8+4))
!            write(6,550) RpRs,(bt(j),tflux(j),j=1,nintg)
 550        format(30(F8.5,1X))
         endif

      else
         do j=1,nintg
            tflux(j)=1.0d0
         enddo
      endif

      tmodel=tmodel+sum(tflux)/dble(nintg)-1.0d0

   else
      tmodel=tmodel+0.0d0
   endif
enddo

tmodel=tmodel+1.0d0+zpt !restore zero point

return
end subroutine
