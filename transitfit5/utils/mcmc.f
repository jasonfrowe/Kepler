CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mcmc(npta,aT,aM,aE,aIT,dtype,nfitm,nfit,nplanet,sol,
     .  sol2,serr,err,tmodel,rchi,seed,flag,bchi,ng,dil)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npta,nfit,i,seed,flag,ng,nplanet,j,nfitm,ii
      integer dtype(npta)
      double precision aT(npta),aM(npta),aE(npta),sol(nfitm),aIT(npta),
     .  serr(nfitm,2),tmodel(npta),sol2(nfitm),rchi,ran2,chi1,chi2,
     .  fratio,u,alpha,bchi,gasdev,dil(2),jitter(2),jrvp(2),jrv,rgas,
     .  err(nfitm)
     
      jitter(1)=2.0 !m/s RV stellar jitter
      jitter(2)=1.0 !m/s RV stellar jitter float
C     reformalize jitters for log scale
      jrvp(1)=log(jitter(1))
      jrvp(2)=(log((jitter(1)+jitter(2))/jitter(1)))**2.0d0
      jrv=0.0d0 !default value of jitter.
  
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
c      write(0,*) ng,serr(ng,2)
      if((ng.lt.1).or.(ng.gt.nfit)) goto 15 !make sure we have a valid i
      if(serr(ng,2).eq.0.0d0) goto 15 !if sig is zero we don't change
      
      j=ng-10*((ng-9)/10)
      write(0,*) ng
      if(ng.eq.1)then
        sol2(ng)=gasdev(seed)*serr(ng,2)+sol(ng)
        do 22 i=1,nplanet
            ii=11+10*(i-1)
 23         sol2(ii)=sol2(ng)*serr(ii,1)+err(ii)+gasdev(seed)*serr(ii,2)
            if((sol2(ii).lt.0.0d0).or.(sol2(ii).gt.1.0d0)) then
                write(0,*) ng,ii,sol2(ii)
c                goto 23
                flag=1
                return
            endif
            ii=12+10*(i-1)
 24         sol2(ii)=sol2(ng)*serr(ii,1)+err(ii)+gasdev(seed)*serr(ii,2)
            if(sol2(ii).lt.0.0d0) then
                write(0,*) ng,ii,sol2(ii)
c                goto 24
                flag=1
                return
            endif
 22     continue
      elseif((j.eq.9).and.(ng.gt.8))then
        sol2(ng)=gasdev(seed)*serr(ng,2)+sol(ng)
        sol2(ng+1)=gasdev(seed)*serr(ng+1,2)+sol(ng+1)
  
      elseif((j.eq.10).and.(ng.gt.8))then
        sol2(ng-1)=gasdev(seed)*serr(ng-1,2)+sol(ng-1)
        sol2(ng)=gasdev(seed)*serr(ng,2)+sol(ng)
c        write(0,*) 'ng: ',ng,ng-1
        
      elseif((j.eq.11).and.(ng.gt.8))then
        sol2(1)=gasdev(seed)*serr(1,2)+sol(1)
        do 25 i=1,nplanet
            ii=11+10*(i-1)
 26         sol2(ii)=sol2(1)*serr(ii,1)+err(ii)+gasdev(seed)*serr(ii,2)
            if((sol2(ii).lt.0.0d0).or.(sol2(ii).gt.1.0d0)) then
                write(0,*) ng,ii,sol2(ii)
c                goto 26
                flag=1
                return
            endif
            ii=12+10*(i-1)
 27         sol2(ii)=sol2(1)*serr(ii,1)+err(ii)+gasdev(seed)*serr(ii,2)
            if(sol2(ii).lt.0.0d0) then
                write(0,*) ng,ii,sol2(ii)
c                goto 27
                flag=1
                return
            endif
 25     continue
c 18     sol2(ng)=gasdev(seed)*serr(ng,2)+sol(ng)
c        if(sol2(ng).lt.0.0d0) goto 18
c        if(sol2(ng).gt.1.0d0) goto 18
c 19     sol2(ng+1)=gasdev(seed)*serr(ng+1,2)+sol(ng+1)
c        if(sol2(ng+1).lt.0.0d0) goto 19
cc        write(0,*) 'ng: ',ng,ng+1
        
      elseif((j.eq.12).and.(ng.gt.8))then
        sol2(1)=gasdev(seed)*serr(1,2)+sol(1)
        do 28 i=1,nplanet
            ii=11+10*(i-1)
 29         sol2(ii)=sol2(1)*serr(ii,1)+err(ii)+gasdev(seed)*serr(ii,2)
            if((sol2(ii).lt.0.0d0).or.(sol2(ii).gt.1.0d0)) then
                write(0,*) ng,ii,sol2(ii)
c                goto 29
                flag=1
                return
            endif
            ii=12+10*(i-1)
 30         sol2(ii)=sol2(1)*serr(ii,1)+err(ii)+gasdev(seed)*serr(ii,2)
            if(sol2(ii).lt.0.0d0) then
                write(0,*) ng,ii,sol2(ii)
c                goto 30
                flag=1
                return
            endif
 28     continue
c 20     sol2(ng-1)=gasdev(seed)*serr(ng-1,2)+sol(ng-1)
c        if(sol2(ng-1).lt.0.0d0) goto 20
c        if(sol2(ng-1).gt.1.0d0) goto 20
c 21     sol2(ng)=gasdev(seed)*serr(ng,2)+sol(ng)
c        if(sol2(ng).lt.0.0d0) goto 21
cc        write(0,*) 'ng: ',ng,ng-1
        
      else
        sol2(ng)=gasdev(seed)*serr(ng,2)+sol(ng)
      endif
c      if(sol2(3).lt.0.0d0) then
c        flag=1
c        return
c      endif
C     Shuffle dilution with Gaussian for every chain element
 17   sol2(6)=gasdev(seed)*dil(2)+dil(1)
      if(sol2(6).lt.0.0d0) goto 17
      if(jitter(2).gt.0.0d0)then
        rgas=gasdev(seed)
        jrv=exp(rgas*jrvp(2)+jrvp(1))
c        write(0,*) "jrv",jrv,rgas*jitter(2)+jitter(1)
c        read(5,*)
      endif
      
C     Get chi-squared for both models
      call transitmodel(nfit,nplanet,sol,npta,aT,aIT,tmodel,dtype)
      chi1=0.0d0
      do 11 i=1,npta
        if(dtype(i).eq.1)then
            chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i)+
     .          jrv*jrv)
        else
            chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
        endif
c        chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
c        chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/tmodel(i)
 11   continue
      chi1=chi1*bchi
 
      call transitmodel(nfit,nplanet,sol2,npta,aT,aIT,tmodel,dtype)
      chi2=0.0d0
      do 12 i=1,npta
        if(dtype(i).eq.1)then
            chi2=chi2+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i)+
     .          jrv*jrv)
        else
            chi2=chi2+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
        endif
c        chi2=chi2+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
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