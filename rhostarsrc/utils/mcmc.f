      subroutine mcmc(seed,nsol,sol,sol2,dsol,nmodelmax,tage,tTeff,
     .  tlogL,trad,y2,Teff,Tefferr,Z,Zerr,L,Lerr,logg,loggerr,nfrho,
     .  rhoi,rhoierr,rhoin,yprho,flag,ng,rchi,nfitm,gscale,ngcor,
     .  ngcorsub,ngprob,ngprobsub,nbuffer,nb,nmov,buffer,nup,
     .  npars,nupdate,corscale,nacor,nacorsub,naprob,naprobsub,tmcore,
     .  aotu,tdloggdt)
C     Implementation of efficient and automatic deMCMC
C     Jason Rowe - jasonfrowe@gmail.com
      implicit none
      integer seed,nmodelmax,nmodel,flag,ng,i,ntry,ntrymax,nfrho,nsol,
     .  nfit
      double precision sol(nsol),dsol(3),tage(nmodelmax),
     .  tTeff(nmodelmax),tlogL(nmodelmax),trad(nmodelmax),gasdev,
     .  sol2(nsol),yp1,ypn,y2(nmodelmax),Teffn2,radn2,fratio,chi1,chi2,
     .  u,alpha,rchi,ran2,Teff,Tefferr,rhoi,rhoierr(9),rhoin(9),
     .  yprho(9),Rsun,Msun,Ms,Rs,Teffn1,radn1,rhon1,rhon2,Pi,drho,dsig,
     .  trho(nmodelmax),tdrhodt(nmodelmax),drhodt1,drhodt2,Z,Zerr(2),
     .  L,Lerr,logL1,logL2,logg,loggerr,logg1,logg2,tmcore(nmodelmax),
     .  mcore,mcore2,mcoresig,tagemax,tagemin,tdloggdt(nmodelmax),
     .  dloggdt1,dloggdt2
      integer nfitm,ngcor(nfitm),ngcorsub(nfitm),ngprob(nfitm),
     .  ngprobsub(nfitm),nbuffer,nupcor,nb,nmov,nup,npars,nupdate,
     .  nsel,nsel2,nacor,nacorsub,naprob,naprobsub,ageprior
      double precision gscale(nfitm),gcorsub,mcmctype,
     .  buffer(nfitm,nbuffer),corscale,corsub,aotu
      character*80 cout
      
      ageprior=0 ! ageprior=0 use prior, =1 do not use prior
      
      flag=0 !flag=0 means accept, flag=1 means reject
      nupcor=0 !default we do not update Gibbs samplers or Correlations
      
      Pi=acos(-1.d0)!define Pi and 2*Pi
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Msun=1.9891d30 !kg  mass of Sun

c      write(0,*) "gt1",sol(1),sol(3)
      call gettrack(sol(1),sol(3),nmodelmax,nmodel,tage,tTeff,tlogL,
     .  trad,trho,tdrhodt,tmcore,tdloggdt,Tefferr,loggerr)
c      write(0,*) "gt2"
     
      yp1=1.0e30
      ypn=1.0e30
c      call spline(tage,tTeff,nmodel,yp1,ypn,y2)
c      call splint(tage,tTeff,y2,nmodel,sol(2),Teffn1)
      call lininterp(tage,tTeff,nmodel,sol(2),Teffn1)
c      call spline(tage,trad,nmodel,yp1,ypn,y2)
c      call splint(tage,trad,y2,nmodel,sol(2),radn1)
      call lininterp(tage,trad,nmodel,sol(2),radn1) 
c      call spline(tage,tdrhodt,nmodel,yp1,ypn,y2)
c      call splint(tage,tdrhodt,y2,nmodel,sol(2),drhodt1)
      call lininterp(tage,tdrhodt,nmodel,sol(2),drhodt1) 
      call lininterp(tage,tdloggdt,nmodel,sol(2),dloggdt1)      
      call lininterp(tage,tlogL,nmodel,sol(2),logL1)
      call lininterp(tage,tmcore,nmodel,sol(2),mcore)
      
      sol(4)=radn1
      sol(6)=Teffn1
      sol(7)=logL1
      Ms=sol(1)*Msun
      Rs=radn1*Rsun
      
      chi1=0.0d0 !initialize chi-square stat to zero
      rhon1=Ms/(4.0d0/3.0d0*pi*Rs*Rs*Rs)
      logg1=log10(6.67259E-8*sol(1)*1.989E33/
     .  (6.9599E10*6.9599E10*radn1*radn1))
      sol(5)=rhon1
c        write(6,*) Teffn1,radn1,rhon1
      if(nfrho.eq.0)then  !if nfrho=0, then we are fitting rho_*
        drho=rhon1-rhoi
c        call splint(rhoierr,rhoin,yprho,9,drho,dsig)
        call getrhosig(rhoierr,rhoin,9,drho,dsig)
c        write(6,*) "dsig:",dsig,(Teffn1-Teff)/Tefferr
        chi1=chi1+dsig*dsig
      endif
      if(Tefferr.gt.0.0d0)then !if Tefferr>0, then we are fiting Teff
        chi1=chi1+((Teffn1-Teff)/Tefferr)**2.0d0
      endif
      if(Zerr(1).gt.0.0d0)then !if Zerr>0, then we are fitting Z
c        write(0,*) Z,sol(3),sol(3)-Z
        if(sol(3)-Z.ge.0.0d0)then
            chi1=chi1+((sol(3)-Z)/Zerr(1))**2.0d0
        else
            chi1=chi1+((sol(3)-Z)/Zerr(2))**2.0d0
        endif
      endif
      if(Lerr.gt.0.0d0)then !if Lerr>0, then we are fitting L
        chi1=chi1+((10**logL1-L)/Lerr)**2.0d0
      endif
      if(loggerr.gt.0.0d0)then!If log(g)_err > 0, then we are fit log(g)
        chi1=chi1+((logg1-logg)/loggerr)**2.0d0
      endif
c      if(ageprior.eq.0) chi1=chi1+mcoresig(tage,tmcore,nmodel,sol(2),
c     .  mcore)

      ntrymax=30
      ntry=0
      nfit=3
      do 14 i=1,nfit
        sol2(i)=sol(i)
 14   continue
      if(dsol(3).eq.0.0d0) nfit=2

      mcmctype=ran2(seed)

      if((nmov.lt.nbuffer).or.(mcmctype.le.0.5))then            
        ng=int(ran2(seed)*dble(nfit)+1.0d0)
 15     continue
        if(ng.eq.1)then
 11         sol2(1)=gasdev(seed)*dsol(1)*gscale(1)+sol(1) !Mass
c            if(sol2(1).lt.0.4) goto 11
c            if(sol2(1).gt.5.2) goto 11
        elseif(ng.eq.2)then
 12         sol2(2)=gasdev(seed)*dsol(2)*gscale(2)+sol(2) !Age
c            if(sol2(2).lt.0.0) goto 12
c            if(sol2(2).gt.tage(nmodel)) goto 12
        elseif(ng.eq.3)then
 13         sol2(3)=gasdev(seed)*dsol(3)*gscale(3)+sol(3) !Z
c            if(sol2(3).lt.0.00001) goto 13
c            if(sol2(3).gt.0.08) goto 13
        endif
      else
 21     nsel=int(ran2(seed)*dble(nbuffer-1)+1.0d0)
        nsel2=int(ran2(seed)*dble(nbuffer-1)+1.0d0)
        do 20 i=1,nfit
            sol2(i)=sol(i)+(buffer(i,nsel2)-buffer(i,nsel))*corscale
 20     continue

C       THIS IS NOT TECHNICALLY CORRECT, WE SHOULD EXIT WITH FLAG=1 IF
C       WE FALL OUTSIDE BOUNDS.
CR        if(sol2(1).lt.0.4) goto 21 !Keep Mass inside isochrone bounds
CR        if(sol2(1).gt.5.2) goto 21
CR        if(sol2(2).lt.0.0) goto 21 !Keep Age inside bounds
CR        if(sol2(2).gt.tage(nmodel)) goto 21
CR        if(sol2(3).lt.0.00001) goto 21 !Keep Z inside bounds
CR        if(sol2(3).gt.0.08) goto 21        
c        if((sol2(1).lt.0.4).or.(sol2(1).gt.5.2).or.(sol2(2).lt.0.0).or.
c     .      (sol2(2).gt.80.0d0).or.(sol2(3).lt.0.00001).or.
c     .      (sol2(3).gt.0.08))then
c                flag=1!jump past getting new parameters and reject chain
c        endif
c        write(0,*) (sol2(i),i=1,nfit)
      endif
c        ntry=ntry+1

      if((sol2(1).lt.0.4).or.(sol2(1).gt.5.2).or.(sol2(2).lt.0.0).or.
     .  (sol2(2).gt.aotu).or.(sol2(3).lt.0.00001).or.
     .  (sol2(3).gt.0.08))then
            flag=1!jump past getting new parameters and reject chain
      endif

      
      if(flag.eq.0)then
c      write(0,*) "gt3",sol(1),sol(3)
c      write(0,*) mcmctype,corscale
        call gettrack(sol2(1),sol2(3),nmodelmax,nmodel,tage,tTeff,tlogL,
     .      trad,trho,tdrhodt,tmcore,tdloggdt,Tefferr,loggerr)
c      write(0,*) "gt4"
     
        yp1=1.0e30
        ypn=1.0e30
c      call spline(tage,tTeff,nmodel,yp1,ypn,y2)
c      call splint(tage,tTeff,y2,nmodel,sol2(2),Teffn2)
        call lininterp(tage,tTeff,nmodel,sol2(2),Teffn2)
CR      if(Teffn2.le.0.0d0)then !make sure Teff is sane
CR        if((nmov.lt.nbuffer).or.(mcmctype.le.0.5))then
CR            goto 15
CR        else
CR           goto 21
CR        endif
CR     endif

c      call spline(tage,trad,nmodel,yp1,ypn,y2)
c      call splint(tage,trad,y2,nmodel,sol2(2),radn2)
        call lininterp(tage,trad,nmodel,sol2(2),radn2)
CR      if(radn2.le.0.0d0) then !make sure radius is sane
CR        if((nmov.lt.nbuffer).or.(mcmctype.le.0.5))then
CR            goto 15
CR        else
CR            goto 21
CR        endif
CR      endif

c      call spline(tage,tdrhodt,nmodel,yp1,ypn,y2)
c      call splint(tage,tdrhodt,y2,nmodel,sol2(2),drhodt2)
        call lininterp(tage,tdrhodt,nmodel,sol2(2),drhodt2)
        call lininterp(tage,tdloggdt,nmodel,sol2(2),dloggdt2)  
        call lininterp(tage,tlogL,nmodel,sol2(2),logL2)
        call lininterp(tage,tmcore,nmodel,sol2(2),mcore2)
      
        sol2(4)=radn2
        sol2(6)=Teffn2
        sol2(7)=logL2
        Ms=sol2(1)*Msun
        Rs=radn2*Rsun
      
        chi2=0.0d0 !initialize chi-square stat to zero
        rhon2=Ms/(4.0d0/3.0d0*pi*Rs*Rs*Rs)
        logg2=log10(6.67259E-8*sol2(1)*1.989E33/
     .      (6.9599E10*6.9599E10*radn2*radn2))
C     Make sure mean density is between 0 and 20 g/cm^3
CR      if((rhon2.le.0.0).or.(rhon2.gt.20000.0)) goto 15
        sol2(5)=rhon2
        if(nfrho.eq.0)then  !if nfrho=0, then we are fitting rho_*
c       write(6,*) Teffn2,radn2,rhon2
            drho=rhon2-rhoi
c       call splint(rhoierr,rhoin,yprho,9,drho,dsig)
            call getrhosig(rhoierr,rhoin,9,drho,dsig)
c       write(6,*) "dsig:",dsig,(Teffn2-Teff)/Tefferr
            chi2=chi2+dsig*dsig
        endif
        if(Tefferr.gt.0.0d0)then !if Tefferr>0, then we are fiting Teff
            chi2=chi2+((Teffn2-Teff)/Tefferr)**2.0d0
        endif
        if(Zerr(1).gt.0.0d0)then !if Zerr>0, then we are fitting Z
            if(sol2(3)-Z.ge.0.0d0)then
                chi2=chi2+((sol2(3)-Z)/Zerr(1))**2.0d0
            else
                chi2=chi2+((sol2(3)-Z)/Zerr(2))**2.0d0
            endif
        endif
c        if(Zerr.gt.0.0d0)then !if Zerr>0, then we are fitting Z
c            chi2=chi2+((sol2(3)-Z)/Zerr)**2.0d0
c        endif
        if(Lerr.gt.0.0d0)then !if Lerr>0, then we are fitting L
            chi2=chi2+((10**logL2-L)/Lerr)**2.0d0
        endif
        if(loggerr.gt.0.0d0)then!If log(g)_err>0, then we are fit log(g)
            chi2=chi2+((logg2-logg)/loggerr)**2.0d0
        endif
c        if(ageprior.eq.0) chi2=chi2+mcoresig(tage,tmcore,nmodel,
c     .      sol2(2),mcore2)
     
c        write(0,*) sol2(1),sol2(2),
c     .      mcoresig(tage,tmcore,nmodel,sol2(2),mcore2)
      
c      write(0,*) logg1,logg,loggerr
c      read(5,*)
      
c      chi2=dsig*dsig+((Teffn2-Teff)/Tefferr)**2.0d0+
c     .  ((sol2(3)-Z)/Zerr)**2.0d0
c      write(0,*) sol(2),sol2(2),abs(drhodt1/drhodt2)
      
c        write(0,*) sol(2),sol2(2)
c        write(0,*) 1.0d0-mcore,1.0d0-mcore2
c        write(0,*) (1.0d0-mcore)/(1.0d0-mcore2)
c        read(5,*)
c        fratio=exp(0.5d0*(chi1-chi2))!/abs((1.0d0-mcore)/(1.0d0-mcore2))
C                                               !/abs(drhodt1/drhodt2)
c        write(0,*) dloggdt1,dloggdt2,abs(dloggdt1/dloggdt2)
c        write(0,*) sol(2),sol2(2)
        if(ageprior.eq.0)then
            fratio=exp(0.5d0*(chi1-chi2))/abs(dloggdt1/dloggdt2)
        else
            fratio=exp(0.5d0*(chi1-chi2))
        endif
        u=ran2(seed)
        alpha=min(fratio,1.0d0)
        if(u.le.alpha)then
            flag=0
            rchi=chi2  !accept jump
        else
            flag=1 
            rchi=chi1  !reject jump
        endif

      endif

C     Move all the GOTO 15 statements here. (CR comments above)
C     reject jump if radius is below zero
      if(sol2(4).le.0.0) flag=1
C     reject jump if density is below zero.
      if((sol2(5).le.0.0).or.(sol2(5).gt.20000.0)) flag=1 
C     reject jump if Teff is below zero
      if(sol2(6).le.0.0) flag=1
      
      tagemax=tage(1)
      tagemin=tage(1)
      do 22 i=2,nmodel
        tagemax=max(tage(i),tagemax)
        tagemin=min(tage(i),tagemin)
 22   continue
      if((sol2(2).gt.tagemax).or.(sol2(2).lt.tagemin)) then
        flag=1 !if age is outside range, reject!
c        write(0,*) "age:",sol2(1),sol2(2),tagemax
      endif


C     Fill Buffer (initial fill)
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
c            npars=npars+2  !keep a larger history..
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
      
c      if((flag.eq.1).and.(ntry.lt.ntrymax))goto 15

c      nupcor=0
c      if(flag.eq.0)then
c        nb=nb+1  !count number of accepted trials
cc        write(0,*) nb
c        if(nb.eq.nbuffer)then
c            nupcor=1
c            nb=1
c        endif
c      endif

      if((nmov.lt.nbuffer).or.(mcmctype.le.0.5))then
        if(flag.eq.0)then
            ngcor(ng)=ngcor(ng)+1
            ngcorsub(ng)=ngcorsub(ng)+1
        endif
        ngprob(ng)=ngprob(ng)+1
        ngprobsub(ng)=ngprobsub(ng)+1
      else
        if(flag.eq.0)then
            nacor=nacor+1
            nacorsub=nacorsub+1
        endif
        naprob=naprob+1
        naprobsub=naprobsub+1
      endif
      
      if(nupcor.eq.1)then
        do 16 i=1,nfitm
            if(ngprob(i).ne.ngprobsub(i))then
                gcorsub=dble(ngcor(i)-ngcorsub(i))/
     .              dble(ngprob(i)-ngprobsub(i))
                gscale(i)=gscale(i)*(0.75*(gcorsub+0.01)/(0.25*
     .              (1.0d0-gcorsub+0.01)))**0.25d0
            endif
            ngcorsub(i)=0 !reset counter
            ngprobsub(i)=0 !reset counter
 16     continue
        if(naprob.ne.naprobsub)then
            corsub=dble(nacor-nacorsub)/
     .          dble(naprob-naprobsub)
            corscale=corscale*(0.75*(corsub+0.01)/(0.25*
     .          (1.0d0-corsub+0.01)))**0.25d0
c            write(0,*) naprob,naprobsub,corsub
c            write(0,*)  nacor,nacorsub
        endif
        nacorsub=0 !reset counter
        naprobsub=0 !reset counter
        write(0,*)
        write(0,*) "corscale"
        write(0,500) corscale 
        write(0,*) "gscale"
        write(0,500) (gscale(i),i=1,nfitm)
        write(0,*) 
 500    format(20(F5.1,1X)) 
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
      integer len,iwhich
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
