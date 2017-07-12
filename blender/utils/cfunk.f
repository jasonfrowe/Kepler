CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function cfunk(x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,i,nfit,npars
      parameter(nmax=600000,nfit=16)
      integer ipars(nfit)
      double precision x(nfit),time(nmax),mag(nmax),merr(nmax),
     .  itime(nmax),sol(nfit),sol2(nfit),lml,Pi,r2pi,tmodel(nmax),chi,
     .  serr(nfit,2),gauss,t1,t2,t3
      common /Fitting/ npars,ipars,npt,time,mag,merr,itime,sol,serr,
     .  sol2,tmodel

      Pi=acos(-1.0)
      r2pi=Sqrt(2.0*Pi)

      do 10 i=1,nfit
         sol2(i)=sol(i)
 10   continue
      do 13 i=1,npars
         sol2(ipars(i))=x(i)
 13   continue
 501  format(13(1PE10.3,1X))
c      call transitmodel(npt,time,tmodel,sol2)
      call transitmodel(npt,time,itime,tmodel,nfit,sol2)
      
      lml=0.
      chi=0.
      do 11 i=1,npt
         lml=lml+log(1.0/(merr(i)*r2pi))-
     .      (mag(i)-tmodel(i))*(mag(i)-tmodel(i))/(2.0*merr(i)*merr(i))
c         chi=chi+(mag(i)-tmodel(i))*(mag(i)-tmodel(i))/(merr(i)*merr(i))
 11   continue
      do 12 i=1,nfit
         if(serr(i,2).gt.0)then
            t1=serr(i,2)
            t2=serr(i,1)
            t3=sol2(i)
            gauss=exp(-(t3-t2)*(t3-t2)/(2.0d0*t1*t1))
            lml=lml+log(gauss)
c            lml=lml+log(rgauss(serr(i,2),serr(i,1),sol2(i)))
         endif
 12   continue

c      cfunk=1.0/lml
      cfunk=-lml
c      write(6,*) cfunk,x
c      read(5,*)
c      cfunk=chi

      return
      end
 
      
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      double precision function rgauss(sig,mu,x)
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      implicit none
c      double precision sig,mu,x,temp(2),Pi
c      Pi=acos(-1.d0)
c      
cc      temp(1)=sig*sqrt(2.0*Pi)
c      temp(1)=1.0
c      temp(2)=(x-mu)*(x-mu)/(2.0d0*sig*sig)
cc      write(6,*) "t",temp(2)
c      rgauss=exp(-temp(2))/temp(1)
c      
c      return
c      end