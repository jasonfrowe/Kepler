CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dmcmc(npta,aT,aM,aE,aIT,dtype,nfitm,nfit,nplanet,sol,
     .  sol2,serr,err,tmodel,rchi,seed,flag,bchi,ng,dil,nup,nb,nmov,
     .  nbuffer,buffer,npars,corscale,nacor,nacorsub,naprob,naprobsub,
     .  gscale,ngcor,ngcorsub,ngprob,ngprobsub,nupdate)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npta,nfitm,nfit,nplanet,seed,flag,ng,i,nbuffer,nmov,nup,
     .  nb,ng2,npars,nsel,nsel2,j,nupcor,nupdate
      integer dtype(npta),nacor(nfitm,nfitm),nacorsub(nfitm,nfitm),
     .  ngcor(nfitm),ngcorsub(nfitm),naprob(nfitm,nfitm),
     .  naprobsub(nfitm,nfitm),ngprob(nfitm),ngprobsub(nfitm)
      double precision aT(npta),aM(npta),aE(npta),sol(nfitm),aIT(npta),
     .  sol2(nfitm),serr(nfitm,2),err(nfitm),tmodel(npta),rchi,bchi,
     .  dil(2),jitter(2),jrvp(2),jrv,ran2,dm,gasdev,rgas,chi1,chi2,
     .  fratio,u,alpha,buffer(nfitm,nbuffer),mcmctype,dif1,dif2,
     .  corscale(nfitm,nfitm),acorsub,gscale(nfitm),gcorsub
      character*80 cout 
     
      nupcor=0 !default we do not update Gibbs samplers or Correlations
     
      jitter(1)=2.0 !m/s RV stellar jitter
      jitter(2)=1.0 !m/s RV stellar jitter float
C     reformalize jitters for log scale
      jrvp(1)=log(jitter(1))
      jrvp(2)=(log((jitter(1)+jitter(2))/jitter(1)))**2.0d0
      jrv=0.0d0 !default value of jitter.
      
      do 16 i=1,nfit
        sol2(i)=sol(i)
  16  continue
  
  15  ng=int(ran2(seed)*dble(nfit-1)+1.0d0)
c      write(0,*) ng,serr(ng,2)
      if((ng.lt.1).or.(ng.gt.nfit)) goto 15 !make sure we have a valid i
      if(ng.eq.6) goto 15 !dilution is not fitted.
      if(serr(ng,2).eq.0.0d0) goto 15 !if sig is zero we don't change
      
      mcmctype=ran2(seed)
      if((nmov.lt.nbuffer).or.(mcmctype.le.0.5))then
        sol2(ng)=gasdev(seed)*serr(ng,2)*gscale(ng)+sol(ng)
      else
 20     ng2=int(ran2(seed)*dble(nfit-1)+1.0d0)
        if(ng2.eq.ng) goto 20 !make sure we have 2 different ngs 
        if((ng2.lt.1).or.(ng2.gt.nfit)) goto 20 !check for valid i
        if(ng2.eq.6) goto 20 !dilution is not fitted.
        if(serr(ng2,2).eq.0.0d0) goto 20 !if sig is zero we don't change
        nsel=int(ran2(seed)*dble(nbuffer-1)+1.0d0)
        nsel2=int(ran2(seed)*dble(nbuffer-1)+1.0d0)
        dif1=buffer(ng,nsel2)-buffer(ng,nsel)
        dif2=buffer(ng2,nsel2)-buffer(ng2,nsel)
        if((dif1.eq.0.0d0).and.(dif2.eq.0.0d0)) goto 20 !make sol move
        sol2(ng)=sol(ng)+dif1*corscale(ng,ng2)
        sol2(ng2)=sol(ng2)+dif2*corscale(ng,ng2)
      endif
      
 17   sol2(6)=gasdev(seed)*dil(2)+dil(1)  !update dilution.
      if(sol2(6).lt.0.0d0) goto 17  
      if(jitter(2).gt.0.0d0)then   !add RV jitter
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
        flag=0  !accept chain
        rchi=chi2
      else
        flag=1  !reject chain
        rchi=chi1
      endif
c      write(0,*) floor(mcmctype+0.5),flag

C     Filling buffer      
      if((flag.eq.0).and.(nmov.lt.nbuffer))then
        nmov=nmov+1
        do 18 i=1,nfit
            buffer(i,nmov)=sol2(i)
 18     continue
        write(cout,501) "nmov:",nmov
 501    format(A5,I4)
        call ovrwrt(cout,2)
      endif

C     Update buffer if accepted
      if((flag.eq.0).and.(nup.ge.npars).and.(nmov.eq.nbuffer).and.
     .   (mcmctype.le.0.5))then
        nup=0
        nb=nb+1
        if(nb.eq.nbuffer)then    
            nupcor=1
            nupdate=nupdate+1
            npars=npars+2  !keep a larger history..
        endif
        if(nb.gt.nbuffer) nb=1
        write(cout,502) "nb: ",nb,nupdate
 502    format(A4,I4,1X,I4)
        call ovrwrt(cout,2)
        do 19 i=1,nfit
            buffer(i,nb)=sol2(i)
 19     continue
      elseif((flag.eq.0).and.(nup.lt.npars).and.(nmov.eq.nbuffer).and.
     .   (mcmctype.le.0.5))then
        nup=nup+1
      endif
          
C     scale counts
      if((nmov.lt.nbuffer).or.(mcmctype.le.0.5))then
        if(flag.eq.0)then
            ngcor(ng)=ngcor(ng)+1
            if(nupdate.gt.0) ngcorsub(ng)=ngcorsub(ng)+1
        endif
        ngprob(ng)=ngprob(ng)+1
        if(nupdate.gt.0) ngprobsub(ng)=ngprobsub(ng)+1
      else
        if(flag.eq.0)then
            nacor(ng,ng2)=nacor(ng,ng2)+1
            if(nupdate.gt.0) then
                nacorsub(ng,ng2)=nacorsub(ng,ng2)+1
c                write(0,*) ng,ng2,nacorsub(ng,ng2)
            endif
        endif
        naprob(ng,ng2)=naprob(ng,ng2)+1
        if(nupdate.gt.0) then
            naprobsub(ng,ng2)=naprobsub(ng,ng2)+1
c            write(0,*) ng,ng2,naprobsub(ng,ng2)
        endif
      endif
     
C     Update corscale and gscale
      if((nupcor.eq.1).and.(nupdate.gt.1))then
        write(0,*) "Updating corscale:",nupdate-1
        do 21 i=1,18!nfitm
            do 22 j=1,18!nfitm
                if(naprob(i,j).ne.naprobsub(i,j))then
                    acorsub=dble(nacor(i,j)-nacorsub(i,j))/
     .                  dble(naprob(i,j)-naprobsub(i,j))
                    corscale(i,j)=corscale(i,j)*(0.75*(acorsub+0.01)/
     .                  (0.25*(1.0d0-acorsub+0.01)))**0.25d0
                endif
                nacorsub(i,j)=0 !reset counter
                naprobsub(i,j)=0 !reset counter
 22         continue
            if(ngprob(i).ne.ngprobsub(i))then
                gcorsub=dble(ngcor(i)-ngcorsub(i))/
     .              dble(ngprob(i)-ngprobsub(i))
                gscale(i)=gscale(i)*(0.75*(gcorsub+0.01)/(0.25*
     .              (1.0d0-gcorsub+0.01)))**0.25d0
            endif
            ngcorsub(i)=0 !reset counter
            ngprobsub(i)=0 !reset counter
 21     continue
        do 23 i=1,18
            write(0,500) (corscale(i,j),j=1,18)
 500        format(20(F5.2,1X))
 23     continue
        write(0,*) "gscale"
        write(0,500) (gscale(i),i=1,18)
        write(0,*)
        nupcor=0 !done updating
      endif
 
      return
      end