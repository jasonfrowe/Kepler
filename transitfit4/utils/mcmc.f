CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mcmc(npta,aT,aM,aE,aIT,dtype,nfit,sol,sol2,serr,tmodel,
     .  rchi,seed,flag,bchi,ng,dil)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npta,nfit,i,seed,flag,ng
      integer dtype(npta)
      double precision aT(npta),aM(npta),aE(npta),sol(nfit),aIT(npta),
     .  serr(nfit,2),tmodel(npta),sol2(nfit),rchi,ran2,chi1,chi2,
     .  fratio,u,alpha,bchi,gasdev,dil(2)
  
c 14   continue    
c      do 10 i=1,nfit
c 13      continue
cc         sol2(i)=serr(i,1)+serr(i,2)*(2.0d0*ran2(seed)-1.0d0)
cc         if(i.eq.3)then
cc            sol2(i)=serr(i,1)+serr(i,2)*(2.0d0*ran2(seed)-1.0d0)
cc         else
cc            sol2(i)=gasdev(seed)*serr(i,2)+serr(i,1)
cc         endif
c         sol2(i)=gasdev(seed)*serr(i,2)+sol(i)
c         if(i.eq.3)then
c            if(sol2(i).lt.0.0d0) goto 13
c            if(sol2(i).gt.1.0d0) goto 13
c         endif
c         if((i.eq.7).or.(i.eq.8))then
c            if(sol2(i).lt.-0.95d0) goto 13
c            if(sol2(i).gt. 0.95d0) goto 13
c         endif
c 10   continue
c      if(sqrt(sol(7)*sol(7)+sol(8)*sol(8)).gt.0.95d0) goto 14

      do 16 i=1,nfit
        sol2(i)=sol(i)
  16  continue
  
  15  ng=int(ran2(seed)*dble(nfit-1)+1.0d0)
      if((ng.lt.1).or.(ng.gt.nfit-1)) goto 15 !make sure we have a valid i
      if(serr(ng,2).eq.0.0d0) goto 15 !if sig is zero we don't change
      if((ng.ge.1).and.(ng.le.2))then
        sol2(1)=gasdev(seed)*serr(1,2)+sol(1)
        sol2(2)=gasdev(seed)*serr(2,2)+sol(2)
      elseif((ng.eq.3).or.(ng.eq.4))then
 18     sol2(3)=gasdev(seed)*serr(3,2)+sol(3)
        if(sol2(3).lt.0.0d0) goto 18
        if(sol2(3).gt.1.0d0) goto 18
 19     sol2(4)=gasdev(seed)*serr(4,2)+sol(4)
c 19     sol2(4)=gasdev(seed)*serr(4,2)+0.1168+sol(3)*0.015
        if(sol2(4).lt.0.0d0) goto 19
      elseif((ng.eq.5).or.(ng.eq.8))then
        sol2(5)=gasdev(seed)*serr(5,2)+sol(5)
        sol2(8)=gasdev(seed)*serr(8,2)+sol(8)
      elseif(ng.eq.6)then
        sol2(6)=gasdev(seed)*serr(6,2)+sol(6)
      elseif((ng.eq.7).or.(ng.eq.10))then
        sol2(7)=gasdev(seed)*serr(7,2)+sol(7)
        sol2(10)=gasdev(seed)*serr(10,2)+sol(10)
      elseif(ng.eq.9)then
        sol2(9)=gasdev(seed)*serr(9,2)+sol(9)
      else
        sol2(ng)=gasdev(seed)*serr(ng,2)+sol(ng)
      endif
      if(sol2(3).lt.0.0d0) then
        flag=1
        return
      endif
C     Shuffle dilution with Gaussian for every chain element
 17   sol2(18)=gasdev(seed)*dil(2)+dil(1)
      if(sol2(18).lt.0.0d0) goto 17
      
C     Get chi-squared for both models
      call transitmodel(nfit,sol,npta,aT,aIT,tmodel,dtype)
      chi1=0.0d0
      do 11 i=1,npta
        chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
c        chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/tmodel(i)
 11   continue
      chi1=chi1*bchi
 
      call transitmodel(nfit,sol2,npta,aT,aIT,tmodel,dtype)
      chi2=0.0d0
      do 12 i=1,npta
        chi2=chi2+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
c        chi2=chi2+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/tmodel(i)
c        write(0,*) aM(i),tmodel(i)
 12   continue
      chi2=chi2*bchi
 
      fratio=exp(0.5d0*(chi1-chi2))
      u=ran2(seed)
      alpha=min(fratio,1.0d0)
      if(u.le.alpha)then
        flag=0
        rchi=chi2
      else
        flag=1
        rchi=chi1
      endif
c      write(0,*) flag,fratio,u
c      write(0,*) alpha,chi1,chi2
c      write(0,*) sol2(7),sol2(8)
c      read(5,*)
     
      return
      end