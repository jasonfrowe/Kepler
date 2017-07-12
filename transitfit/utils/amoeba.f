CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER iter,mp,ndim,np,NMAX,ITMAX
      REAL*8 ftol,p(mp,np),y(mp),funk,TINY
      PARAMETER (NMAX=20,ITMAX=500,TINY=1.d-10)
      character*80 tmpc
C     Maximum allowed dimensions and function evaluations, and a small 
C     number.
      EXTERNAL funk
C     USES amotry,funk
C     Multidimensional minimization of the function funk(x) where 
C     x(1:ndim) is a vector in ndim dimensions, by the downhill simplex 
C     method of Nelder and Mead. The matrix p(1:ndim+1,1:ndim) is input.
C     Its ndim+1 rows are ndim-dimensional vectors which are the 
C     vertices of the starting simplex. Also input is the vector 
C     y(1:ndim+1), whose components must be pre-initialized to the 
C     values of funk evaluated at the ndim+1 vertices (rows) of p; and 
C     ftol the fractional convergence tolerance to be achieved in the 
C     function value (n.b.!). On output, p and y will have been reset to
C     ndim+1 new points all within ftol of a minimum function value, and 
C     iter gives the number of function evaluations taken.
      INTEGER i,ihi,ilo,inhi,j,m,n
      REAL*8 rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
      iter=0

C     Enter here when starting or have just overall contracted.
 1    do 12 n=1,ndim
         sum=0. !Recompute psum.
         do 11 m=1,ndim+1
            sum=sum+p(m,n)
 11      continue
         psum(n)=sum
 12   continue
 2    ilo=1 !Enter here when have just changed a single point.
      if (y(1).gt.y(2)) then !Determine which point is the highest 
C                       (worst), next-highest,and lowest (best)
         ihi=1 
         inhi=2
      else
         ihi=2
         inhi=1
      endif
      do 13 i=1,ndim+1 !by looping over the points in the simplex.
         if(y(i).le.y(ilo)) ilo=i
         if(y(i).gt.y(ihi)) then
            inhi=ihi
            ihi=i
         else if(y(i).gt.y(inhi)) then
            if(i.ne.ihi) inhi=i
         endif
 13   continue
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
C     Compute the fractional range from highest to lowest and return if 
C     satisfactory.
      if (rtol.lt.ftol) then 
C     If returning, put best point and value in slot 1.
         swap=y(1)
         y(1)=y(ilo)
         y(ilo)=swap
         do 14 n=1,ndim
            swap=p(1,n)
            p(1,n)=p(ilo,n)
            p(ilo,n)=swap
 14      continue
         return
      endif
      if (iter.ge.ITMAX) return!pause "ITMAX exceeded in amoeba"
      iter=iter+2
C     Begin a new iteration. First extrapolate by a factor −1 through 
C     the face of the simplex across from the high point, i.e., reflect 
C     the simplex from the high point.
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0)
      if (ytry.le.y(ilo)) then
C     Gives a result better than the best point, so try an additional 
C     extrapolation by a factor 2.
         ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0)
      else if (ytry.ge.y(inhi)) then
C     The reflected point is worse than the second-highest, so look for 
C     an intermediate lower point, i.e., do a one-dimensional 
C     contraction.
         ysave=y(ihi)
         ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5)
         if (ytry.ge.ysave) then 
C     Can’t seem to get rid of that high point. Better contract
C     around the lowest (best) point.
            do 16 i=1,ndim+1
               if(i.ne.ilo)then
                  do 15 j=1,ndim
                     psum(j)=0.5*(p(i,j)+p(ilo,j))
                     p(i,j)=psum(j)
 15               continue
                  y(i)=funk(psum)
               endif
 16         continue
            iter=iter+ndim !Keep track of function evaluations.
            if(iter.gt.495) write(0,*) ytry,ysave,iter
            goto 1 
C     Go back for the test of doneness and the next iteration.
         endif
      else
         iter=iter-1 !Correct the evaluation count.
      endif
      if(mod(iter,10).eq.0)then
        write(tmpc,500) "Iteration # ", iter
 500    format(A12,I4)
        call ovrwrt(tmpc,2)
      endif
      goto 2
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ihi,mp,ndim,np,NMAX
      REAL*8 fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=20)
      EXTERNAL funk
C     USES funk
C     Extrapolates by a factor fac through the face of the simplex 
C     across from the high point, tries it, and replaces the high point 
C     if the new point is better.
      INTEGER j
      REAL*8 fac1,fac2,ytry,ptry(NMAX)
      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
         ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
 11   continue
      ytry=funk(ptry) !Evaluate the function at the trial point.
      if (ytry.lt.y(ihi)) then 
C     If it’s better than the highest, then replace the highest.
         y(ihi)=ytry
         do 12 j=1,ndim
            psum(j)=psum(j)-p(ihi,j)+ptry(j)
            p(ihi,j)=ptry(j)
 12      continue
      endif
      amotry=ytry
      return
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ovrwrt (line, iwhich)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Taken from the DAOPhot code by Stetson.
C     This is a cool routine for writing lots of info to the screen
C     and keeping everything on the same line.
C
      character*(*) line
      character*79 output
      integer len,iwhich
      if (iwhich .eq. 1) then
         write (0,1) line
    1    format (a)
      else if (iwhich .eq. 2) then
         if (len(line) .lt. 79) then
            output = ' '
            output = line
            write (0,2) output, char(13), char(13)
            write (0,2) output, char(13), char(13)
            write (0,2) output, char(13), char(13)
    2       format (a, 2a1, $)
         else
            write (0,2) line, char(13), char(13)
         end if
      else if (iwhich .eq. 3) then
         write (0,3) line
    3    format (a)
      else
         write (0,4) line, char(13), char(13)
    4    format (/a, 2a1, $)
         write (0,2) line, char(13), char(13)
      end if
      return
      end