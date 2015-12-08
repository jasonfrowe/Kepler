CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine shaperm(npt,nfit,time,mag,merr,res,pers,avg,a,panx,
     .     pany,plot,w,nstar,Keplertime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Given fourier coefficents, this shape is removed from the data.
      implicit none
      integer npt,nfit,na,nstarmax,panx,pany,plot
      parameter(nstarmax=300)
      double precision time(npt),mag(npt),res(npt),pers(nstarmax),avg,
     .  a(nfit),merr(npt),Keplertime
      integer i,j,k,nstar,i2,cc
      double precision pi,w(nstarmax),temp,arg
      character*80 title

      pi=3.141592654
      do 5 i=1,nstar
         w(i)=(2.0*pi)/pers(i)
 5    continue
      na=nstar*nfit*2+nstar

      j=(na-nstar)/nstar
      do 10 i2=1,npt
         temp=avg
         do 20 k=1,nstar
            cc=(j+1)*(k-1)+2
            do 21 i=cc,cc+j-2,2
               arg=dble((i-cc+2)/2)*a(cc-1)*time(i2)+a(i+1)
               temp=temp+a(i)*cos(arg)
 21         continue
 20      continue
         res(i2)=mag(i2)-temp
 10   continue
      
      title="Residuals"
      if(plot.eq.1) call plotdata(npt,time,res,merr,-1.0d0,-1.0d0,panx,
     .  pany,Keplertime,title)
      
      return
      end