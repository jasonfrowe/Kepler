      subroutine getbit(qflag,qfbit)
      implicit none
      integer qflag,qfbit(14),i
      character*14 cbin

      write(cbin,500) qflat
 500  format(B14.1)
      do 10 i=1,14
         qfbit(i)=0
         read(cbin(i:i),*,end=10) qfbit(i)
 10   continue

      write(6,501) (qfbit(i),i=1,14)
 501  format(14I1)

      return
      end
