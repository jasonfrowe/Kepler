subroutine octiming(nbodies,nintg,xpos,ypos,zpos,sol,opos,tc,tcalc,told,&
 bold,ntmid,tmid)
use precision
implicit none
integer, dimension(:) :: ntmid
real(double), dimension(:,:) :: tmid
integer nbodies,i,j,nintg,np,k
integer, dimension(:) :: tc
real(double) :: rstar,RpRs,LU2,b,rt1,rt2,s1,s2,b1,b2,bmin,fourpid3,&
  rhostar,Mstar
real(double), dimension(:,:) :: xpos,ypos,zpos,told,bold
real(double), dimension(:) :: sol,opos,tcalc
!vars for output of timing measurements
integer :: itprint

itprint=0 

LU2=LU*LU
!rstar=(sol(12)/(4.0d0/3.0d0*Pi*sol(1)*1000.0d0)*Mearth)**(1.0d0/3.0d0)
fourpid3=4.0d0*Pi/3.0d0
rhostar=abs(sol(1)*1000.0) !mean stellar density (kg/m^3)
Mstar=abs(sol(12)) !mass of central object (MEarth)
rstar=(Mearth*Mstar/(fourpid3*rhostar))**(1.0d0/3.0d0)

do i=2,nbodies
   np=7+7*(i-1)
   bmin=abs(sol(np+3)) !minimum impact parameter. 
   RpRs=abs(sol(np+4))
   do j=1,nintg
      !check if we are now on the transit side
      if((ypos(i,j)-ypos(1,j).gt.0.0d0).and.(opos(i)-ypos(1,j).le.0.0d0))then 
         tc(i)=0
      endif
      if(ypos(i,j)-ypos(1,j).le.0.0d0)then
         tc(i)=-1 !on eclipse side.. might as well skip and wait...
      endif

      if(tc(i).ge.0)then  !if on side of transit then check for transit
!     see if we have a transit
         b=(zpos(i,j)-zpos(1,j))*(zpos(i,j)-zpos(1,j))+ &
            (xpos(i,j)-xpos(1,j))*(xpos(i,j)-xpos(1,j))
         b=sqrt(LU2*b)/rstar
         if(b.lt.1.0d0+RpRs)then !if we have a transit, then store time
            tc(i)=tc(i)+1
            told(i,tc(i))=tcalc(j) !time in current transit
            bold(i,tc(i))=b !projection of radii
         endif

         if(tc(i).eq.5)then !we need 5 transit points to do the calculation
            rt1=(bold(i,4)-bold(i,5))/(bold(i,2)-bold(i,3)) !compare slopes
            rt2=(bold(i,3)-bold(i,4))/(bold(i,1)-bold(i,2)) !compare slopes
!            write(0,*) bold(1)-bold(2),bold(2)-bold(3),bold(3)-bold(4), &
!             bold(4)-bold(5)
!            if(i.eq.2) write(0,*) tcalc(j),rt1,rt2
            if((rt1.lt.0).or.(rt2.lt.0))then !see if we passed the min.
!              write(0,*) "o-c!"
               tc(i)=-1 ! no point in searching for more minimum
               s1=(bold(i,1)-bold(i,2))/(told(i,1)-told(i,2)) !calculate slopes
               s2=(bold(i,4)-bold(i,5))/(told(i,4)-told(i,5))
               
               b1=bold(i,1)-s1*told(i,1)
               b2=bold(i,5)-s2*told(i,5)
               ntmid(i)=ntmid(i)+1
               tmid(i,ntmid(i))=(b2-b1)/(s1-s2) !estimate of mid transit time.
!               if(i.eq.2) write(0,*) "tmid:",ntmid(i),tmid(i,ntmid(i))
               
               if(itprint.eq.1)then
                  if (i.eq.4) then
                  !   write(0,501) told(i,1),bold(i,1),b1,bmin
                  !   write(0,501) told(i,2),bold(i,2)
                  !   write(0,501) told(i,3),bold(i,3)
                  !   write(0,501) told(i,4),bold(i,4)
                  !   write(0,501) told(i,5),bold(i,5),b2,bmin
                     if (ntmid(i)-1 > 0) then
                        write(0,501) tmid(i,ntmid(i)),tmid(i,ntmid(i))-tmid(i,ntmid(i)-1)
                     endif
                  !   read(5,* )
                  endif
               endif
               501  format(78(1PE13.6,1X))

               
            else
               tc(i)=4 ! no minimum.. so delete earliest point
               do k=1,4
                  told(i,k)=told(i,k+1)
                  bold(i,k)=bold(i,k+1)
               enddo
            endif
         endif

      endif

      opos(i)=ypos(i,j)
   enddo

   !write(0,501) "percal",i,ntmid(i)
   !501 format(A6,1X,I3,1X,I3,1X,28(1PE15.8,1X))

enddo


return
end subroutine octiming
