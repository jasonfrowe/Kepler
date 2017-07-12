!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine f(neq,t,y,ydot)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use precision
      implicit none
      integer neq,npt,i,j
      double precision t,y(neq),ydot(neq),tmp(3),R!,Rmin

      npt=nbod
!      write(0,*) "npt:",npt
      do 10 i=1,npt
!         write(0,*) "G:",G
!         write(0,*) "M",i,m(i)
!     Okay.  for x,y,z the derivatives are simply vx,vy,vz.  simple.
         ydot(6*i-5)=y(6*i-2)
         ydot(6*i-4)=y(6*i-1)
         ydot(6*i-3)=y(6*i)
!     Initialize tmp array
         do 12 j=1,3
            tmp(j)=0.0d0
 12      continue
         do 11 j=1,npt
            if((j.ne.i).and.(m(j).gt.0.0d0))then
               R=(y(6*i-5)-y(6*j-5))*(y(6*i-5)-y(6*j-5))+  &
                 (y(6*i-4)-y(6*j-4))*(y(6*i-4)-y(6*j-4))+  &
                 (y(6*i-3)-y(6*j-3))*(y(6*i-3)-y(6*j-3))
!               Rmin=min(R,Rmin)
               R=R**(3.0d0/2.0d0)
               tmp(1)=tmp(1)+G*m(j)*(y(6*j-5)-y(6*i-5))/R
               tmp(2)=tmp(2)+G*m(j)*(y(6*j-4)-y(6*i-4))/R
               tmp(3)=tmp(3)+G*m(j)*(y(6*j-3)-y(6*i-3))/R
            endif
 11      continue
         ydot(6*i-2)=tmp(1)
         ydot(6*i-1)=tmp(2)
         ydot(6*i)=tmp(3)
!         write(6,*)(ydot(j),j=6*i-5,6*i)
 10   continue

      return
      end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE JAC(VEC,MAXNEQ)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
      INTEGER MAXNEQ
      double precision VEC(maxneq)

      RETURN
      END
