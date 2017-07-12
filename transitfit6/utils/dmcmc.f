CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dmcmc(npta,aT,aM,aE,aIT,dtype,nfitm,nfit,nplanet,sol,
     .  sol2,serr,err,tmodel,rchi,seed,flag,bchi,ng,dil,nup,nb,nmov,
     .  nbuffer,buffer,npars,corscale,nacor,nacorsub,naprob,naprobsub,
     .  gscale,ngcor,ngcorsub,ngprob,ngprobsub,nupdate,gratio,nfrho,
     .  rhoi,rhoierr,rhoin,nrad,ntheta,rbb)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npta,nfitm,nfit,nplanet,seed,flag,ng,i,nbuffer,nmov,nup,
     .  nb,ng2,npars,nsel,nsel2,j,nupcor,nupdate,nplot,nrad,ntheta
      integer dtype(npta),nacor,nacorsub,
     .  ngcor(nfitm),ngcorsub(nfitm),naprob,
     .  naprobsub,ngprob(nfitm),ngprobsub(nfitm),njit
      real rbb(3,4)
      double precision aT(npta),aM(npta),aE(npta),sol(nfitm),aIT(npta),
     .  sol2(nfitm),serr(nfitm,2),err(nfitm),tmodel(npta),rchi,bchi,
     .  dil(2),jitter(2),jrvp(2),jrv,ran2,dm,gasdev,rgas,chi1,chi2,
     .  fratio,u,alpha,buffer(nfitm,nbuffer),mcmctype,dif1,dif2,
     .  corscale,acorsub,gscale(nfitm),gcorsub,gratio(nfitm),echeck,b
      integer nfrho
      double precision rhoi,rhoierr(9),rhoin(9),dsig,drho
      character*80 cout 
      
      nplot=0 !nplot=0, plot some stuff, nplot=1, do not plot
     
      nupcor=0 !default we do not update Gibbs samplers or Correlations
     
      njit=0 !if njit=0, no RV jitter. 
      jitter(1)=2.0 !m/s RV stellar jitter
      jitter(2)=1.0 !m/s RV stellar jitter float
C     reformalize jitters for log scale
      jrvp(1)=log(jitter(1))
      jrvp(2)=(log((jitter(1)+jitter(2))/jitter(1)))**2.0d0
      jrv=0.0d0 !default value of jitter.
      
      do 16 i=1,nfit
        sol2(i)=sol(i)
  16  continue
  
  15  ng=int(ran2(seed)*dble(nfit)+1.0d0)
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
      if(sol2(6).lt.0.0d0) goto 17  
      if(jitter(2).gt.0.0d0)then   !add RV jitter
        rgas=gasdev(seed)
        if(njit.eq.0)then
            jrv=0.0d0
        else
            jrv=exp(rgas*jrvp(2)+jrvp(1))
        endif
      endif
      
      if(rchi.lt.0.0)then !see if this is the first run)
C       Get chi-squared for both models
        call transitmodel(nrad,ntheta,nfit,nplanet,sol,npta,aT,aIT,
     .  dtype,tmodel,rbb,nplot)
        chi1=0.0d0
        do 11 i=1,npta
            if(dtype(i).eq.1)then
                chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*
     .            aE(i)+jrv*jrv)
            else
                chi1=chi1+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*
     .            aE(i))
c               write(0,*) i,chi1,tmodel(i)
            endif
 11     continue
        chi1=chi1*bchi      
        if(nfrho.eq.0)then  !if nfrho=0, then we are fitting rho_*
            drho=1.0d3*sol(1)-rhoi
            call getrhosig(rhoierr,rhoin,9,drho,dsig)
            chi1=chi1+dsig*dsig
        endif
      else
        chi1=rchi
      endif
 
c      call transitmodel(nfit,nplanet,sol2,npta,aT,aIT,tmodel,dtype)
      call transitmodel(nrad,ntheta,nfit,nplanet,sol2,npta,aT,aIT,dtype,
     .  tmodel,rbb,nplot)
      chi2=0.0d0
      do 12 i=1,npta
        if(dtype(i).eq.1)then
            chi2=chi2+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i)+
     .          jrv*jrv)
        else
            chi2=chi2+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
c            write(0,*) i,chi2,tmodel(i)
        endif
c        chi2=chi2+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
c        chi2=chi2+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/tmodel(i)
c        write(0,*) aM(i),tmodel(i)
 12   continue
      chi2=chi2*bchi
 
      if(nfrho.eq.0)then  !if nfrho=0, then we are fitting rho_*
        drho=1.0d3*sol2(1)-rhoi
c        call splint(rhoierr,rhoin,yprho,9,drho,dsig)
        call getrhosig(rhoierr,rhoin,9,drho,dsig)
c        write(6,*) "dsig:",dsig,(Teffn1-Teff)/Tefferr
        chi2=chi2+dsig*dsig
c        write(0,*) 1.0d3*sol2(1),rhoi,dsig
c        read(5,*)
      endif 
 
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
      
      if((sol2(1).lt.0.001).or.(sol2(1).gt.8.0)) flag=1 !density
      if((sol2(12).lt.0.001).or.(sol2(12).ge.1.0)) flag=1 !f
      do 25 i=1,nplanet
        b=sol2(18+11*(i-1))
        if((b.lt.0.0).or.
     .     (b.ge.1.0+sol2(19+11*(i-1)))) flag=1 !impact parameter
        if((sol2(19+11*(i-1)).lt.0.0).or.
     .     (sol2(19+11*(i-1)).gt.1.0)) flag=1 !r/R*
        if((sol2(20+11*(i-1)).lt.-1.0).or.
     .     (sol2(20+11*(i-1)).gt.1.0)) flag=1 !ecosw
        if((sol2(21+11*(i-1)).lt.-1.0).or.
     .     (sol2(21+11*(i-1)).gt.1.0)) flag=1 !esinw
        echeck=sol2(20+11*(i-1))**2 + sol2(21+11*(i-1))**2
        if((echeck.lt.0.0).or.(echeck.ge.1.0)) flag=1 !e
        if(sol2(22+11*(i-1)).lt.0.0) flag=1  !K
 25   continue
      
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
cc        write(0,*) "Updating corscale:",nupdate-1
cc        if(naprob.ne.naprobsub)then
cc            acorsub=dble(nacor-nacorsub)/
cc     .          dble(naprob-naprobsub)
cc            corscale=corscale*(0.75*(acorsub+0.01)/(0.25*
cc     .          (1.0d0-acorsub+0.01)))**0.25d0
cc        endif
        write(0,*) "corscale"
        write(0,*) nacorsub,naprobsub,dble(nacorsub)/dble(naprobsub)
        nacorsub=0 !reset counter
        naprobsub=0 !reset counter

        do 21 i=1,18!nfitm
            if(ngprob(i).ne.ngprobsub(i))then
                gcorsub=dble(ngcor(i)-ngcorsub(i))/
     .              dble(ngprob(i)-ngprobsub(i))
                gscale(i)=gscale(i)*(0.75*(gcorsub+0.01)/(0.25*
     .              (1.0d0-gcorsub+0.01)))**0.25d0
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
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ovrwrt (line, iwhich)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Taken from the DAOPhot code by Stetson.
C     This is a cool routine for writing lots of info to the screen
C     and keeping everything on the same line.
C
      character*(*) line
      character*79 output
      integer len
      if (iwhich .eq. 1) then
         write (0,1) line
    1    format (a)
      else if (iwhich .eq. 2) then
         if (len(line) .lt. 79) then
            output = ' '
            output = line
            write (0,2) output, char(13), char(13)
            write (0,2) output, char(13), char(13)
            write (0,2) output, char(13), char(13)
    2       format (a, 2a1, $)
         else
            write (0,2) line, char(13), char(13)
         end if
      else if (iwhich .eq. 3) then
         write (0,3) line
    3    format (a)
      else
         write (0,4) line, char(13), char(13)
    4    format (/a, 2a1, $)
         write (0,2) line, char(13), char(13)
      end if
      return
      end