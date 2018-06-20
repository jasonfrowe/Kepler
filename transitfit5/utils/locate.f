      subroutine locate(array,length,x,idx)
      integer idx,length
      real*8 x,array(length)
c     Given an array and a value x, returns idx
c     such that x is between array(idx) and array(idx+1).

      integer low,mid,up
      low=0                             !initialize limits
      up=length+1
 10   if(up-low.gt.1)then
         mid=(up+low)/2                 !compute midpoint
         if((array(length).ge.array(1)).eqv.(x.ge.array(mid)))then
            low=mid             !replace either lower or upper limit
         else
            up=mid
         endif
         goto 10
      endif
      if(x.eq.array(1))then
         idx=1
      else if(x.eq.array(length))then
         idx=length-1
      else
         idx=low
      endif
      return
      end
