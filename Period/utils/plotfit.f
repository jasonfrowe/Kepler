CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotfit(npt,nfit,time,mag,pers,avg,a,w,nstar)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nfit,nstarmax,nstar,i,j,k,i2,cc,nmax,na,steps
      parameter(nstarmax=300,nmax=2000000,steps=2000)
      real px(nmax),py(nmax)
      double precision time(npt),mag(npt),pers(nstarmax),avg,a(nfit),
     .  w(nstarmax),pi,temp,arg,tt,maxt,mint

      pi=3.141592654
      do 5 i=1,nstar
         w(i)=(2.0*pi)/pers(i)
 5    continue
      na=nstar*nfit*2+nstar

      maxt=-99.9e30
      mint= 99.9e30
      do 7 i=1,npt
         maxt=max(maxt,time(i))
         mint=min(mint,time(i))
 7    continue

      j=(na-nstar)/nstar
      do 10 i2=1,steps
         tt=real(i2)*(maxt-mint)/real(steps)+mint
         temp=avg
c         write(6,*) avg
         do 20 k=1,nstar
            cc=(j+1)*(k-1)+2
            do 21 i=cc,cc+j-2,2
c               write(6,*) cc-1,(i-cc+2)/2,a(cc-1),a(i),a(i+1)
               arg=real((i-cc+2)/2)*a(cc-1)*tt+a(i+1)
               temp=temp+a(i)*cos(arg)
 21         continue
 20      continue
         px(i2)=real(tt)
         py(i2)=real(temp)
 10   continue
      call pgsci(2)
      call pgline(steps,px,py)
      call pgsci(1)

      return
      end
