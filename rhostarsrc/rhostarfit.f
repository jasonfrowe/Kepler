      program rhostarfit
c     find best fit values for Teff, logg and [Fe/H]
      implicit none 
      integer iargc,m,n,nmax,info,nfit,lwa,nmass,nmodelmax,nmodel
      parameter(nmax=3,nfit=3)
      parameter(nmass=36,nmodelmax=1000)
      parameter(lwa=nmax*nfit+5*nmax*nfit)
      integer iwa(nmax),ia(nfit)
      double precision Teff,logg,FeH,Tefferr,loggerr,feherr,Z,Zerr,dum,
     .  aotu,mass,rad,age,x(nfit),eps,fvec(nmax),wa(lwa),tol,radn1,
     .  tage(nmodelmax),tTeff(nmodelmax),tlogL(nmodelmax),
     .  trad(nmodelmax),trho(nmodelmax),tdrhodt(nmodelmax),
     .  tmcore(nmodelmax),tdloggdt(nmodelmax),Teffn1,logL1,rhon1,logg1,
     .  Pi,Rsun,Msun,Ms,Rs,covar(nfit,nfit),alpha(nfit,nfit),a(nfit),
     .  y(nfit),sig(nfit),dchistop
      character*80 cline
      common /Fitting/ Teff,Tefferr,logg,loggerr,Z,Zerr,aotu
c      common /Fitting2/ Tefferr,loggerr
      
      aotu=14.0d0 !Approximate Age Of The Universe.
      eps=0.0001 !small error -a way of handling error=0 for fixed value
      Pi=acos(-1.d0)!define Pi and 2*Pi
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Msun=1.9891d30 !kg  mass of Sun

C     Read in Teff and error      
      if(iargc().lt.6) goto 901
      call getarg(1,cline)
      read(cline,*) Teff
      if((Teff.lt.3000.0).or.(Teff.gt.100000)) goto 902
      call getarg(2,cline)
      read(cline,*) Tefferr
      if(Tefferr.lt.0.0) goto 905
      if(Tefferr.eq.0.0) Tefferr=eps
      
C     Read in logg and error
      call getarg(3,cline)
      read(cline,*) logg
      if((logg.le.0.0).or.(logg.gt.10.0)) goto 903
      call getarg(4,cline)
      read(cline,*) loggerr
      if(loggerr.lt.0.0) goto 905
      if(loggerr.eq.0.0) loggerr=eps

C     Read in metallicity and error      
      call getarg(5,cline)
      read(cline,*) FeH
      if((FeH.le.-5.0).or.(FeH.gt.2.0)) goto 904
      call getarg(6,cline)
      read(cline,*) feherr
      if(feherr.lt.0.0) goto 905
      if(feherr.eq.0.0) feherr=eps

C     Convert [Fe/H] to Z
      call feh2z(FeH,Z)
      dum=FeH+FeHerr
      call feh2z(dum,Zerr)
      Zerr=Zerr-Z

C     Get initial guess for mass and radius      
      call massscan(mass,rad,age,Teff,Tefferr,logg,loggerr,Z,aotu)
c      write(6,500) mass,rad,age
 
C     minimize chi-square
c      call fitmassrad(nmax,nfit,mass,age,lwa,iwa,fvec,wa,x)

      dchistop=1.0d-10
      call fitmassrad2(nmax,nfit,x,y,sig,a,ia,covar,alpha,Teff,
     .  Tefferr,logg,loggerr,FeH,FeHerr,mass,rad,Z,dchistop)



C     Now lets get the radius based on mass,age and Z      
      call gettrack(x(1),x(3),nmodelmax,nmodel,tage,tTeff,tlogL,trad,
     .      trho,tdrhodt,tmcore,tdloggdt,Tefferr,loggerr)
      call lininterp(tage,trad,nmodel,x(2),radn1) !radius
      call lininterp(tage,tTeff,nmodel,x(2),Teffn1) !get Teff
      call lininterp(tage,tlogL,nmodel,x(2),logL1) !get logL
      
      Ms=x(1)*Msun
      Rs=radn1*Rsun 
      rhon1=Ms/(4.0d0/3.0d0*pi*Rs*Rs*Rs)/1.0d3 !mean stellar density
      logg1=log10(6.67259E-8*x(1)*1.989E33/
     .  (6.9599E10*6.9599E10*radn1*radn1))  !logg (cgs)
      
C                  mass,radius,Z,Teff,logL,Age,logg,rhostar
      write(6,500) x(1),radn1,x(3),Teffn1,logL1,x(2),logg1,rhon1,Teff,
     .  logg,FeH
      
 500  format(F6.3,1X,F7.3,1(1X,F6.3),1X,F6.0,4(1X,F6.3),1X,F6.0,
     .  2(1X,F6.3))
 
      goto 999
 901  write(0,*) "Usage: rhostarfit Teff sig logg sig [Fe/H] sig"
      goto 999
 902  write(0,*) "Teff out of bounds [3000..100000]"
      goto 999
 903  write(0,*) "logg out of bounds [0.0..10.0]"
      goto 999
 904  write(0,*) "[Fe/H] out of bounds [-5.0..2.0]"
      goto 999
 905  write(0,*) "Uncertainties must be >= 0"
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fitmassrad2(nmax,nfit,x,y,sig,a,ia,covar,alpha,Teff,
     .  Tefferr,logg,loggerr,FeH,FeHerr,mass,rad,Z,dchistop)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nmax,nfit,i,ia(nfit),npt,it,itmax
      double precision x(nmax),y(nmax),Teff,Tefferr,logg,loggerr,Z,Zerr,
     .  sig(nmax),FeH,FeHerr,mass,rad,a(nfit),covar(nfit,nfit),
     .  alpha(nfit,nfit),chisq,alamda,bchi,dchi,ochi,alamhi,dchistop
      logical loop
      
      y(1)=Teff
      y(2)=logg
      y(3)=FeH
      sig(1)=Tefferr
      sig(2)=loggerr
      sig(3)=FeHerr
      a(1)=mass
      a(2)=rad
      a(3)=Z
      npt=nmax
      
      do 10 i=1,3
        ia(i)=1.0
 10   continue

      alamhi=100000.
      alamda=-1.0
      loop=.true.
      itmax=50  !takes about 5 iterations to get minimum.
      it=1
      do 30 while(loop)
        call mrqmin(x,y,sig,nmax,a,ia,nfit,covar,alpha,nfit, 
     *     chisq,alamda)
     
         write(0,504) "Best fit for iteration # ",it
 504     format(A25,I4)
         bchi=chisq
         if(it.gt.1)then
            dchi=(ochi-bchi)!/ochi
         else
            dchi=1.0e30
         endif
         write(0,*) "alamda",alamda,dchi
         write(0,503) a(1),a(2),a(3),dchi
         ochi=bchi
 503     format(4(1PE10.3,1X))
         it=it+1
c         if(it.gt.itmax) loop=.false.
c         if((abs(dchi).lt.dchistop).and.(abs(dchi).gt.0.0)) loop=.false.
c         if(abs(dchi).lt.dchistop) then
          if((it.gt.itmax).or.
     .    ((dchi.lt.dchistop).and.(dchi.ge.0.0d0)).or.
     .    (alamda.gt.alamhi))then
            loop=.false.
            alamda=0.0
            call mrqmin(x,y,sig,nmax,a,ia,nfit,covar,alpha,nfit, 
     *          chisq,alamda)
         endif
         write(0,*) "Chi-squared: ",chisq
 30   enddo
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fitmassrad(nmax,nfit,mass,age,lwa,iwa,fvec,wa,x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer m,n,nmax,nfit,lwa,iwa(nmax),info
      double precision x(nfit),mass,age,Teff,Tefferr,logg,loggerr,
     .  Z,Zerr,aotu,tol,fvec(nmax),wa(lwa)
      common /Fitting/ Teff,Tefferr,logg,loggerr,Z,Zerr,aotu
      external fcn
      
C     m is number of functions: Teff,logg,[Fe/H]
      m=nmax
C     n is number of variables: mass, age, Z
      n=nfit
C     x contains initial guess from massscan
      x(1)=mass
      x(2)=age
      x(3)=Z
c     tolerance for fitting
      tol=1.0d-8

C     call least-squares fitter  
      call lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
c      write(6,500) x(1),x(2),x(3)
      write(0,501) "info:",info
 501  format(A5,1X,I2)
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine fcn(m,n,x,fvec,iflag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer m,n,iflag,nmass,nmodelmax,nmodel,i
      parameter(nmass=36,nmodelmax=1000)
      double precision x(n),fvec(m),mass,Ztest,tage(nmodelmax),
     .  tTeff(nmodelmax),tlogL(nmodelmax),trad(nmodelmax),
     .  trho(nmodelmax),tdrhodt(nmodelmax),tmcore(nmodelmax),
     .  tdloggdt(nmodelmax),Teff,Tefferr,logg,loggerr,Z,Zerr,age,
     .  Teffn1,radn1,logg1,aotu,agetest
      common /Fitting/ Teff,Tefferr,logg,loggerr,Z,Zerr,aotu
      
      mass=x(1)
      age=x(2)
      Ztest=x(3)
      
c      write(0,*) x(1),x(2),x(3)
      
      call gettrack(mass,Ztest,nmodelmax,nmodel,tage,tTeff,tlogL,trad,
     .      trho,tdrhodt,tmcore,tdloggdt,Tefferr,loggerr)
     
      agetest=tage(1)
      i=1
      do while(agetest.lt.aotu)
        agetest=tage(i)
        if(agetest.gt.aotu)then
            nmodel=i
        else
            i=i+1
        endif
        if(i.gt.nmodel) agetest=aotu+1.0
      enddo
c      write(0,*) nmodel,tage(nmodel)
c      if(age.gt.tage(nmodel))age=tage(nmodel)-0.01
c      write(0,*) age,aotu,tage(nmodel)
        
     
      call lininterp(tage,tTeff,nmodel,age,Teffn1)
      call lininterp(tage,trad,nmodel,age,radn1) 
      logg1=log10(6.67259E-8*mass*1.989E33/
     .  (6.9599E10*6.9599E10*radn1*radn1))
     
c      write(6,*) mass,radn1,age
 500  format(F8.5,2(1X,F8.5))

      if((age.lt.aotu).and.(age.gt.0.0).and.(mass.gt.0.4).and.
     .  (mass.lt.5.2))then    
        fvec(1)=(Teff-Teffn1)/Tefferr
        fvec(2)=(logg-logg1)/loggerr
        fvec(3)=(Z-Ztest)/Zerr
      else
        fvec(1)=99.9e0
        fvec(2)=99.9e0
        fvec(3)=99.9e0
      endif
      
c      write(0,*) x(1),fvec(1),fvec(2)
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine massscan(massi,radi,agei,Teff,Tefferr,logg,loggerr,
     .  Z,aotu)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,nmass,nmodelmax,nmodel,j,nfrho,k
      parameter(nmass=36,nmodelmax=1000)
      double precision massi,radi,Teff,Tefferr,logg,loggerr,
     .  mass,masses(nmass),tage(nmodelmax),tTeff(nmodelmax),
     .  tlogL(nmodelmax),trad(nmodelmax),trho(nmodelmax),
     .  tdrhodt(nmodelmax),tmcore(nmodelmax),chi1,chiold,rhon1,Teffn1,
     .  logL1,G,Rsun,Msun,Pi,rhocerr,Z,logg1,radn1,dsig,drho,
     .  aotu,agei,tdloggdt(nmodelmax),agetest
      data masses /0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,
     .  1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,
     .  3.6,3.8,4.0,4.1,4.5,5.0,5.2/

      G=6.67259E-11
      Rsun=6.9599E8
      Msun=1.989E30
      Pi=acos(-1.d0)!define Pi and 2*Pi
      

      
      chiold=99.9e30
      do 10 i=1,nmass
        mass=masses(i)
        call gettrack(mass,Z,nmodelmax,nmodel,tage,tTeff,tlogL,trad,
     .      trho,tdrhodt,tmcore,tdloggdt,Tefferr,loggerr)
     

        agetest=tage(1)
        k=1
        do while(agetest.lt.aotu)
            agetest=tage(k)
            if(agetest.gt.aotu)then
                nmodel=k-1
            else
                k=k+1
            endif
            if(k.gt.nmodel) agetest=aotu+1.0
        enddo

     
     
C       We do not scan at boarders, helps chi-square fitter
        do 11 j=2,nmodel-1 
            if(tage(j).gt.aotu) goto 11 !we only ages less that 14Gyr
            rhon1=trho(j)
            Teffn1=tTeff(j)
            logL1=tlogL(j)
            radn1=trad(j)
            logg1=log10(6.67259E-8*mass*1.989E33/
     .          (6.9599E10*6.9599E10*radn1*radn1))
            chi1=0.0d0
            if(Tefferr.gt.0.0d0)then !if Tefferr>0, then fiting Teff
                chi1=chi1+((Teffn1-Teff)/Tefferr)**2.0d0
            endif
            if(loggerr.gt.0.0d0)then!If log(g)_err > 0, then fit log(g)
                chi1=chi1+((logg1-logg)/loggerr)**2.0d0
            endif
            if(chi1.lt.chiold)then
                chiold=chi1
                massi=mass
                radi=radn1
                agei=tage(j)
            endif
 11     continue    
 10   continue
c      write(0,*) chiold,massi,radi

      return
      end     
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine feh2z(FeH,Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (dydz=2.0E0)
      parameter (yp=0.23E0)
      parameter (zp=0.0E0)
      parameter (Zsun=0.0181D0)
      parameter (Xsun=0.7148705D0)
      parameter (FeHa2=-0.217D0)
      parameter (FeHa4=-0.470D0)
      common /head/YYalpfe

      quad(x1,y1,x2,y2,x3,y3,x)=y1*(x2-x)*(x3-x)/((x2-x1)*(x3-x1))
     +                         +y2*(x1-x)*(x3-x)/((x1-x2)*(x3-x2))
     +                         +y3*(x1-x)*(x2-x)/((x1-x3)*(x2-x3))

       FeH0=FeH
     +      -quad(0.0d0,0.0d0,0.3d0,FeHa2,0.6d0,FeHa4,YYalpfe)
C       if(abs(YYalpfe-0.3d0).le.0.1d0)FeH0=FeH-FeHa2
C       if(abs(YYalpfe-0.6d0).le.0.1d0)FeH0=FeH-FeHa4
      ZovX=(10.0d0**FeH0)*Zsun/Xsun
      Z=ZovX*(1.0d0+dydz*zp-yp)/(1.0d0+ZovX*(1.0d0+dydz))
      return
      end
      