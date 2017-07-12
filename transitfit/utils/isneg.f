      !   function isnegative
      !   ----------------
      !   check whether a complex number is a
      !   negative integer
      !
      !   output: 0 if z is zero
      !                   -1 if z is a negative integer
      !                   1  if z is not a negative integer
      !
      integer*4 function isneg(z)
      implicit none
      real*8 az,err
      complex*16 z

      err = 1e-10
      az = cdabs(z)
      if(az.eq.0d0) then
      isneg = 0
        return
      end if
      az = dabs(dreal(z)-dnint(dreal(z)))
      if(az.lt.err.and.dabs(dimag(z)).lt.err.and.dreal(z).lt.0d0) then
        isneg = -1
      else
        isneg = 1
      endif
      return
      end
