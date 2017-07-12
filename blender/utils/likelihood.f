CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function likelihood(npt,tmodel,aT,aM,aE,
     .  aIT,dtype,nfit,sol,serr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Compute likelihood of sol using sol2 and serr as prior 
      implicit none
      integer npt,nfit,i,dtype(npt)
      double precision tmodel(npt),aT(npt),aM(npt),aE(npt),
     .  aIT(npt),sol(nfit),serr(nfit,2),lml,Pi,r2pi,chi,
     .  t1,t2,t3,gauss
      Pi=acos(-1.0)
      r2pi=Sqrt(2.0*Pi)
     
      call transitmodel(npt,aT,aIT,dtype,tmodel,nfit,sol)

      lml=0.
c      write(6,*) "merr:",(merr(i),i=1,npt)
      do 11 i=1,npt
         lml=lml+log(1.0/(aE(i)*r2pi))-
     .      (aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(2.0*aE(i)*aE(i))
c         write(6,*) i,lml,mag(i),tmodel(i)
c         read(5,*)
c         chi=chi+(mag(i)-tmodel(i))*(mag(i)-tmodel(i))/(merr(i)*merr(i))
c         write(6,*) "chi:",chi,mag(i),tmodel(i)
 11   continue
      do 12 i=1,nfit
         if(serr(i,2).gt.0.0d0)then
C           serr(i,2) is the width, serr(i,1) is the center
            t1=serr(i,2)
            t2=serr(i,1)
            t3=sol(i)
            gauss=exp(-(t3-t2)*(t3-t2)/(2.0d0*t1*t1))
            lml=lml+log(gauss)
         endif
 12   continue
 
      likelihood=-lml
c      likelihood=chi
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function chisquared(npt,tmodel,aT,aM,aE,
     .  aIT,dtype,nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,npt,nfit,dtype(npt)
      double precision tmodel(npt),aT(npt),aM(npt),aE(npt),
     .  aIT(npt),sol(nfit),chi
      
      
      call transitmodel(npt,aT,aIT,dtype,tmodel,nfit,sol)
      
      chi=0.0d0
      do 10 i=1,npt
        chi=chi+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
 10   continue
 
      chisquared=chi
      
      return
      end    