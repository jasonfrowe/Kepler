CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dmcmc(npta,aT,aM,aE,aIT,dtype,nfitm,nfit,nplanet,sol,
     .  sol2,serr,err,tmodel,rchi,seed,flag,bchi,ng,dil,nup,nb,nmov,
     .  nbuffer,buffer,npars,corscale,nacor,nacorsub,naprob,naprobsub,
     .  gscale,ngcor,ngcorsub,ngprob,ngprobsub,nupdate,gratio,nfrho,
     .  rhoi,rhoierr,rhoin,nplanetmax,nmax,ntt,tobs,omc,ngs,nas,chiold)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Implementation of deMCMC
C     Jason Rowe - jasonfrowe@gmail.com
      implicit none
      integer npta,nfitm,nfit,nplanet,seed,flag,ng,i,nbuffer,nmov,nup,
     .  nb,ng2,npars,nsel,nsel2,j,nupcor,nupdate
      integer dtype(npta),nacor,nacorsub,
     .  ngcor(nfitm),ngcorsub(nfitm),naprob,
     .  naprobsub,ngprob(nfitm),ngprobsub(nfitm),ngs(nfitm),nas
      double precision aT(npta),aM(npta),aE(npta),sol(nfitm),aIT(npta),
     .  sol2(nfitm),serr(nfitm,2),err(nfitm,2),tmodel(npta),rchi,bchi,
     .  dil(2),jitter(2),jrvp(2),jrv,ran2,dm,gasdev,rgas,chi1,chi2,
     .  fratio,u,alpha,buffer(nfitm,nbuffer),mcmctype,dif1,dif2,
     .  corscale,acorsub,gscale(nfitm),gcorsub,gratio(nfitm),echeck,b,
     .  perprior(2),epoprior(2),chiold,eccn,ecw,esw,w,pi,tpi,rdr
      integer nfrho
      double precision rhoi,rhoierr(9),rhoin(9),dsig,drho,rhostar
      character*80 cout 
      integer nplanetmax,nmax,ntt(nplanetmax)
      double precision tobs(nplanetmax,nmax),omc(nplanetmax,nmax)

      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi

C     Hack for imposing a specific period and epoch
      perprior(1)=3.2357
      perprior(2)=0.0008
      epoprior(1)=3735.17
      epoprior(2)=0.17
     
      nupcor=0 !default we do not update Gibbs samplers or Correlations
     
      jitter(1)=2.0 !m/s RV stellar jitter
      jitter(2)=0.0!1.0 !m/s RV stellar jitter float
C     reformalize jitters for log scale
      jrvp(1)=log(jitter(1))
      jrvp(2)=(log((jitter(1)+jitter(2))/jitter(1)))**2.0d0
      jrv=0.0d0 !default value of jitter.
      
      sol2(1:nfit)=sol(1:nfit)

 15   ng=int(ran2(seed)*dble(nfit)+1.0d0)
c      write(0,*) ng,serr(ng,2)
      if((ng.lt.1).or.(ng.gt.nfit)) goto 15 !make sure we have a valid i
      if(ng.eq.6) goto 15 !dilution is not fitted.
      if(serr(ng,2).eq.0.0d0) goto 15 !if sig is zero we don't change
      
      mcmctype=ran2(seed)
      if((nmov.lt.nbuffer).or.(mcmctype.le.0.5))then
        sol2(ng)=gasdev(seed)*serr(ng,2)*gscale(ng)+sol(ng)
      else
 20     nsel=int(ran2(seed)*dble(nbuffer-1)+1.0d0)
        nsel2=int(ran2(seed)*dble(nbuffer-1)+1.0d0)
        do 24 i=1,nfit
            sol2(i)=sol(i)+(buffer(i,nsel2)-buffer(i,nsel))*corscale
 24     continue
      endif

      
 17   sol2(6)=gasdev(seed)*dil(2)+dil(1)  !update dilution.
c      write(0,*) sol(6),sol2(6),dil(1)
      if(sol2(6).lt.0.0d0) goto 17  
      if(jitter(2).gt.0.0d0)then   !add RV jitter
        rgas=gasdev(seed)
        jrv=exp(rgas*jrvp(2)+jrvp(1))
c        write(0,*) "jrv",jrv,rgas*jitter(2)+jitter(1)
c        read(5,*)
      endif

      if (chiold.le.1.0e-15) then !only compute chi1 for the first iter
C      Get chi-squared for both models
c       call transitmodel(nfit,nplanet,sol,npta,aT,aIT,tmodel,dtype)
       call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npta,aT,aIT,
     .   ntt,tobs,omc,tmodel,dtype)

       chi1=0.0d0
       !$OMP PARALLEL DO REDUCTION(+:chi1)
       do 11 i=1,npta
         if(dtype(i).eq.1)then
             chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i)+
     .           jrv*jrv)
         else
             chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
         endif
c         chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
c         chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/tmodel(i)
 11    continue
       !$OMP END PARALLEL DO
       chi1=chi1*bchi

       if(nfrho.eq.0)then  !if nfrho=0, then we are fitting rho_*
         drho=1.0d3*sol(1)-rhoi
c         call splint(rhoierr,rhoin,yprho,9,drho,dsig)
         call getrhosig(rhoierr,rhoin,9,drho,dsig)
c         write(0,*) "dsig",rhoierr(2),dsig
c         write(6,*) "dsig:",dsig,(Teffn1-Teff)/Tefferr
         chi1=chi1+dsig*dsig
c         write(0,*) 1.0d3*sol(1),rhoi,dsig
       endif

C     Eccentricity constraints...
c      do 26 i=1,8+nplanet*10
c         if(serr(i,1).ne.0)then
c            if(sol(i).gt.serr(i,1))then
c               dsig=(serr(i,1)-sol(i))/err(i,1)
c            else
c               dsig=(serr(i,1)-sol(i))/err(i,2)
c            endif
c            chi1=chi1+dsig*dsig
c         endif
c 26   continue

cC     add in T0 and Period prior (alpha Cen work)
c      dsig=(perprior(1)-sol(10))/perprior(2)
c      chi1=chi1+dsig*dsig
c      dsig=(epoprior(1)-sol(9)-436.0d0*sol(10))/epoprior(2)
c      chi1=chi1+dsig*dsig
cC      write(0,*) dsig,epoprior(1)-sol(9)-436.0d0*sol(10)
      else
         chi1=chiold*bchi
      endif

c      call transitmodel(nfit,nplanet,sol2,npta,aT,aIT,tmodel,dtype)
      call transitmodel(nfit,nplanet,nplanetmax,sol2,nmax,npta,aT,aIT,
     .  ntt,tobs,omc,tmodel,dtype) 
      
      chi2=0.0d0
      !$OMP PARALLEL DO REDUCTION(+:chi2)
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
      !$OMP END PARALLEL DO
      chi2=chi2*bchi

c      chi1=npta !testing  MCMC distributions 
c      chi2=npta
 
      if(nfrho.eq.0)then  !if nfrho=0, then we are fitting rho_*
        drho=1.0d3*sol2(1)-rhoi  !the 10^3 makes sol2(1) kg/m^3
c        call splint(rhoierr,rhoin,yprho,9,drho,dsig)
        call getrhosig(rhoierr,rhoin,9,drho,dsig)
c        write(6,*) "dsig:",dsig,(Teffn1-Teff)/Tefferr
        chi2=chi2+dsig*dsig
c        write(0,*) 1.0d3*sol2(1),rhoi,dsig
c        read(5,*)
      endif 

C     if b>1
      do 28 i=1,nplanet
         b=sol2(11+10*(i-1))
         rdr=sol2(12+10*(i-1))
         if( (b.gt.1.0).and.(rdr.gt.0.0)) then
            chi2=chi2+b*b
         endif
 28   continue

C     prior for Eric Ford
c      do i=1,nplanet
c        b=sol2(11+10*(i-1))
c        rdr=sol2(12+10*(i-1))
c        rhostar=sol2(1)
c        chi2=chi2+rhostar*(1-b*b)/0.13998908451
c      enddo

C     eccentricity constraints
c      do 27 i=1,8+nplanet*10
c         if(serr(i,1).ne.0)then
c            if(sol2(i).gt.serr(i,1))then
c               dsig=(serr(i,1)-sol2(i))/err(i,1)
c            else
c               dsig=(serr(i,1)-sol2(i))/err(i,2)
c            endif
c            chi2=chi2+dsig*dsig
c         endif
c 27   continue

cC     add in T0 and Period prior (alpha Cen work)
c      dsig=(perprior(1)-sol2(10))/perprior(2)
c      chi2=chi2+dsig*dsig
c      dsig=(epoprior(1)-sol2(9)-436.0d0*sol2(10))/epoprior(2)
c      chi2=chi2+dsig*dsig
 
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

c      write(0,*) "rchi:",rchi,npta
      
c      flag=0 !pass everything (for now)

C     bounds for the mean-stellar density
      if((sol2(1).lt.1.0e-5).or.(sol2(1).gt.1000.0)) flag=1 !density

C     bounds for Period (must be positive)
      do 29 i=1,nplanet
         if(sol2(10+10*(i-1)).le.0.0) flag=1 !period 
 29   continue


C     bounds for limb-darkening
      if((sol2(2).eq.0.0).and.(sol2(3).eq.0.0))then  !Kipping Limb-darkening
         if((sol2(4).lt.0.0).or.(sol2(4).gt.1.0).or.
     .    (sol2(5).lt.0.0).or.(sol2(5).gt.1.0))then
            flag=1  !reject model if limb-darkening is out-of-bounds
         endif
      elseif((sol2(4).eq.0.0).and.(sol2(5).eq.0.0))then !Quadratic limb-darkening
         if((sol2(2)+sol2(3).gt.1.0).or.(sol2(2).lt.0).or.
     .    (sol2(2)+2.0d0*sol2(3).lt.0))then
            flag=1  !reject model if limb-darkening is out-of-bounds
         endif
      endif


      ! simple addition to look at MCMC priors (zero-point)
c      if((sol2(8).lt.-5.0e-5).or.(sol2(8).gt.5.0e-5)) flag=1 !zpt

      do 25 i=1,nplanet

      ! simple addition to look at MCMC priors (epoch and period)
c        if((sol2(9+10*(i-1)).lt.113.00).or.
c     .     (sol2(9+10*(i-1)).gt.114.00)) flag=1 !EPO
c        if((sol2(10+10*(i-1)).lt.14.00).or.
c     .     (sol2(10+10*(i-1)).gt.15.00)) flag=1 !Period

        b=sol2(11+10*(i-1)) !changed from b^2 to b
        if((b.lt.0.0).or.
     .     (b.ge.1.0+sol2(12+10*(i-1)))) flag=1 !impact parameter
c     .   (b.gt.1.0)) flag=1
        if((sol2(12+10*(i-1)).lt.0.0).or.
     .     (sol2(12+10*(i-1)).gt.10.0)) flag=1 !r/R*
        if((sol2(13+10*(i-1)).lt.-1.0).or.
     .     (sol2(13+10*(i-1)).gt.1.0)) flag=1 !ecosw
        if((sol2(14+10*(i-1)).lt.-1.0).or.
     .     (sol2(14+10*(i-1)).gt.1.0)) flag=1 !esinw
        echeck=sol2(13+10*(i-1))**2 + sol2(14+10*(i-1))**2
        if((echeck.lt.0.0).or.(echeck.ge.1.0)) flag=1 !e
C        uncomment this line to restrict K > 0
C        if(sol2(15+10*(i-1)).lt.0.0) flag=1  !K
 25   continue
      
c      write(0,*) chi1,chi2,flag
      !save old MCMC value
      if(flag.eq.0)then
         chiold=chi2/bchi
      else
         chiold=chi1/bchi
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
            nacor=nacor+1
            if(nupdate.gt.0) then
                nacorsub=nacorsub+1
c                write(0,*) ng,ng2,nacorsub(ng,ng2)
            endif
        endif
        naprob=naprob+1
        if(nupdate.gt.0) then
            naprobsub=naprobsub+1
c            write(0,*) ng,ng2,naprobsub(ng,ng2)
        endif
      endif
     
C     Update corscale and gscale
      if((nupcor.eq.1).and.(nupdate.gt.1))then
        write(0,*) "Updating corscale:",nupdate-1
c        if(naprob.ne.naprobsub)then
        if((naprob.ne.naprobsub).and.(nas.eq.0))then
        write(0,*) "nas: ",real(nacorsub)/real(naprobsub)
         if((real(nacorsub)/real(naprobsub).lt.0.2).or.(real(nacorsub)/
     .    real(naprobsub).gt.0.3))then
            acorsub=dble(nacor-nacorsub)/
     .          dble(naprob-naprobsub)
            corscale=corscale*(0.75*(acorsub+0.01)/(0.25*
     .          (1.0d0-acorsub+0.01)))**0.25d0
         else
            nas=1 !stop updating scale factor
         endif
        endif
        write(0,*) "corscale"
        write(0,*) nacorsub,naprobsub,dble(nacorsub)/dble(naprobsub)
        nacorsub=0 !reset counter
        naprobsub=0 !reset counter

        do 21 i=1,nfitm
c            if(ngprob(i).ne.ngprobsub(i))then
            if((ngprob(i).ne.ngprobsub(i)).and.(ngs(i).eq.0))then
            write(0,*) "ngs:",i,real(ngcorsub(i))/real(ngprobsub(i))
             if((real(ngcorsub(i))/real(ngprobsub(i)).lt.0.2).or.
     .        (real(ngcorsub(i))/real(ngprobsub(i)).gt.0.3))then
                gcorsub=dble(ngcor(i)-ngcorsub(i))/
     .              dble(ngprob(i)-ngprobsub(i))
                gscale(i)=gscale(i)*(0.75*(gcorsub+0.01)/(0.25*
     .              (1.0d0-gcorsub+0.01)))**0.25d0
             else
               ngs(i)=1 !stop updating scaling factor
             endif
            endif
c            write(0,*) 'ng:',ngcorsub(i),ngprobsub(i)
            gratio(i)=dble(ngcorsub(i))/dble(ngprobsub(i))
            ngcorsub(i)=0 !reset counter
            ngprobsub(i)=0 !reset counter
 21     continue
cc            write(0,500) corscale,acorsub
 500        format(20(F5.2,1X))
        write(0,*) "gscale"
        write(0,500) (gscale(i),i=1,18)
        write(0,500) (gratio(i),i=1,18)
        write(0,*)
        nupcor=0 !done updating
      endif
 
      return
      end
