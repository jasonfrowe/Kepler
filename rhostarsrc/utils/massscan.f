      subroutine massscan(massi,radi,agei,Teff,Tefferr,logg,loggerr,
     .  Psec,L,Lerr,adrsi,adrsierr,rhoc,rhocerr,Z,nfrho,aotu)
      implicit none
      integer i,nmass,nmodelmax,nmodel,j,nfrho
      parameter(nmass=36,nmodelmax=1000)
      double precision massi,radi,Teff,Tefferr,logg,loggerr,Psec,L,Lerr,
     .  adrsi,adrsierr(6),rhoc,rhoi,rhoierr(9),mass,masses(nmass),
     .  tage(nmodelmax),tTeff(nmodelmax),tlogL(nmodelmax),
     .  trad(nmodelmax),trho(nmodelmax),tdrhodt(nmodelmax),
     .  tmcore(nmodelmax),chi1,chiold,rhon1,Teffn1,logL1,G,Rsun,Msun,Pi,
     .  rhocerr,Z,logg1,radn1,dsig,drho,rhoin(9),aotu,agei,
     .  tdloggdt(nmodelmax)
      data masses /0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,
     .  1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,
     .  3.6,3.8,4.0,4.1,4.5,5.0,5.2/
      data rhoin/-1.0d2,-3.0d0,-2.0d0,-1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,
     .  1.0d2/

      G=6.67259E-11
      Rsun=6.9599E8
      Msun=1.989E30
      Pi=acos(-1.d0)!define Pi and 2*Pi
      
C     Calculate rhostar
      rhoi=adrsi**3.0*Pi*3.0d0/(Psec*Psec*G)-rhoc
      rhoierr(1)=-rhoi
      do 15 i=2,4
        rhoierr(i)=(adrsi+adrsierr(i-1))**3.0*Pi*3.0d0/(Psec*Psec*G)
     .      -rhoi-rhoc
        rhoierr(i)=sqrt(rhoierr(i)**2+(dble(5-i)*rhocerr)**2)*
     .      abs(rhoierr(i))/rhoierr(i)        
c        write(0,*) rhoierr(i),adrsierr(i-1)
 15   continue
      rhoierr(5)=0.0d0
      do 16 i=6,8
        rhoierr(i)=(adrsi+adrsierr(i-2))**3.0*Pi*3.0d0/(Psec*Psec*G)
     .      -rhoi-rhoc
        rhoierr(i)=sqrt(rhoierr(i)**2+(dble(i-5)*rhocerr)**2)
c        write(0,*) rhoierr(i),adrsierr(i-2)
 16   continue
      rhoierr(9)=rhoierr(8)*10.0d0!20000.0-rhoi
      
      chiold=99.9e30
      do 10 i=1,nmass
        mass=masses(i)
c        write(0,*) "mass:",mass
        call gettrack(mass,Z,nmodelmax,nmodel,tage,tTeff,tlogL,trad,
     .      trho,tdrhodt,tmcore,tdloggdt,Tefferr,loggerr)
c        write(0,*) "hello?",nmodel
        do 11 j=1,nmodel
c            write(0,*) "hello2",tage(j),aotu
            if(tage(j).gt.aotu) goto 11 !we only ages less that 14Gyr
c            write(0,*) j
            rhon1=trho(j)
            Teffn1=tTeff(j)
            logL1=tlogL(j)
            radn1=trad(j)
            logg1=log10(6.67259E-8*mass*1.989E33/
     .          (6.9599E10*6.9599E10*radn1*radn1))
c            write(0,*) "mr:",mass,radn1
c            write(0,*) "gt:",logg1,Teffn1
            chi1=0.0d0
            if(nfrho.eq.0)then  !if nfrho=0, then we are fitting rho_*
                drho=rhon1-rhoi
c               call splint(rhoierr,rhoin,yprho,9,drho,dsig)
                call getrhosig(rhoierr,rhoin,9,drho,dsig)
c               write(6,*) "dsig:",dsig,(Teffn1-Teff)/Tefferr
                chi1=chi1+dsig*dsig
            endif
c            write(0,*) 1,chi1
            if(Tefferr.gt.0.0d0)then !if Tefferr>0, then fiting Teff
                chi1=chi1+((Teffn1-Teff)/Tefferr)**2.0d0
            endif
c            write(0,*) 2,chi1
            if(Lerr.gt.0.0d0)then !if Lerr>0, then we are fitting L
                chi1=chi1+((10**logL1-L)/Lerr)**2.0d0
            endif
c            write(0,*) 3,chi1
            if(loggerr.gt.0.0d0)then!If log(g)_err > 0, then fit log(g)
                chi1=chi1+((logg1-logg)/loggerr)**2.0d0
            endif
c            write(0,*) 4,chi1
c            if(mass.gt.1.25) write(0,*) chi1,logg1,teffn1
            if(chi1.lt.chiold)then
                chiold=chi1
                massi=mass
                radi=radn1
                agei=tage(j)
c                write(0,*) chi1,massi,radi
c                write(0,*) "Teff",Teffn1,Teff
c                write(0,*) "logg",logg1,logg
            endif
c            read(5,*)
 11     continue    
 10   continue
      write(0,*) chiold,massi,radi

      return
      end