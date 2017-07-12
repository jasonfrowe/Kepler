subroutine octiming(nbodies,nintg,xpos,ypos,zpos,sol,opos,tc,tcalc,told,&
 bold)
use precision
use ocmod
implicit none
integer nbodies,i,j,nintg,np,k
integer, dimension(:) :: tc
real(double) :: rstar,Pi,RpRs,LU2,b,rt1,rt2,m1,m2,b1,b2
real(double), dimension(:,:) :: xpos,ypos,zpos,told,bold
real(double), dimension(:) :: sol,opos,tcalc

Pi=acos(-1.d0)!define Pi
LU2=LU*LU
rstar=(sol(12)/(4.0d0/3.0d0*Pi*sol(1)*1000.0d0)*Mearth)**(1.0d0/3.0d0)

do i=2,nbodies
   np=7+7*(i-1)
   RpRs=abs(sol(np+4))
   do j=1,nintg
      !first we check what side of the star we are on.
      if(ypos(i,j).lt.0.0d0.and.opos(i).ge.0.0)then
         tc(i)=0 !on the transit side.. lets do some computing
      endif
      if(ypos(i,j).ge.0.0d0)then
        tc(i)=-1  !on eclipse side.. might as well skip out
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
               m1=(bold(i,1)-bold(i,2))/(told(i,1)-told(i,2))
               m2=(bold(i,4)-bold(i,5))/(told(i,4)-told(i,5))
               b1=bold(i,1)-m1*told(i,1)
               b2=bold(i,5)-m2*told(i,5)
               ntmid(i)=ntmid(i)+1
               tmid(i,ntmid(i))=(b2-b1)/(m1-m2) !estimate of mid transit time.
               if(i.eq.2) write(0,*) "tmid:",ntmid(i),tmid(i,ntmid(i))
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

enddo


return
end subroutine octiming
