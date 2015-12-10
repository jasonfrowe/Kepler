      double precision function getage(massi,Teff,Tefferr,rhoi,rhoierr,
     .  rhoin,yprho,nmodelmax,nmodel,tTeff,tage,trad,flag,aotu)
      implicit none
      integer nmodelmax,nmodel,i,flag
      double precision Teff,tage(nmodel),tTeff(nmodel),trad(nmodel),
     .  agei,rhoi,rhos,Pi,Rsun,Msun,Ms,Rs,massi,diff,mindiff,minrho,
     .  rhoierr(9),yprho(9),rhoin(9),drho,dsig,Tefferr,aotu
      
      flag=0
      
      Pi=acos(-1.d0)!define Pi and 2*Pi
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Msun=1.9891d30 !kg  mass of Sun
      Ms=massi*Msun
      mindiff=99.9d30
      do 10 i=1,nmodel
        Rs=trad(i)*Rsun
        rhos=Ms/(4.0d0/3.0d0*pi*Rs*Rs*Rs)
        drho=rhos-rhoi
c        call splint(rhoierr,rhoin,yprho,9,drho,dsig)
        call getrhosig(rhoierr,rhoin,9,drho,dsig)
c        diff=sqrt(((Teff-tTeff(i))/Teff)**2 + ((rhoi-rhos)/rhoi)**2)
        diff=dsig*dsig+((tTeff(i)-Teff)/Tefferr)**2.0d0
        if((diff.lt.mindiff).and.(tage(i).lt.aotu))then
            mindiff=diff
            agei=tage(i)
            minrho=rhos
c            write(6,*) mindiff,agei,minrho
c            write(6,*) tTeff(i),trad(i)
        endif
c        write(6,*) diff,tage(i),rhos
c        write(6,*) "T",tTeff(i),trad(i)
 10   continue
      write(0,*) "minrho,agei:",minrho,agei
              
      getage=agei
      
      return
      end
      