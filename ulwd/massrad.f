CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program massradius
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Mass-radius plots for KOI74
      implicit none
      integer nmax,npt,nch,i,nsimp,nc,nfe,nch2,npa,nzp,j
      parameter(nmax=1500)
      real mass(nmax),rad(nmax),masserr(nmax),raderr(nmax),kmass(4),
     .  krad(4),kmasserr(4),kraderr(4),Pi,G,Msun,Mearth,Mjup,Rsun,Rjup,
     .  mch(nmax),rch(nmax),msimp(nmax),rsimp(nmax),dm,mmax,rc(nmax),
     .  mc(nmax),mfe(nmax),rfe(nmax),mch2(nmax),rch2(nmax),pma,pmaerr,
     .  pra,praerr,sma,smaerr,sra,sraerr,teff,palogg(nmax),pamass,
     .  px(nmax),py(nmax),zpm(nmax),zpr(nmax),kmasserr1(4),kmasserr2(4)
      character*80 filename,dumc
      data nch /16/
      data rch /2.273,2.018,1.753,1.558,1.394,1.244,1.094,0.930,0.718,
     .  0.556,0.393,0.299,0.261,0.219,0.189,0.145/
      data mch /0.164,0.224,0.316,0.411,0.512,0.622,0.748,0.899,1.097,
     .  1.235,1.347,1.396,1.411,1.426,1.434,1.444/
      data nc /13/
      data rc /2.111,1.348,1.208,1.067,0.910,0.706,0.548,0.388,0.296,
     .  0.258,0.216,0.208,0.170/
      data mc /0.147,0.488,0.597,0.722,0.872,1.070,1.206,1.318,1.366,
     .  1.381,1.396,1.349,1.174/
      data nfe/21/
      data rfe /2.015,2.102,2.105,2.096,2.004,1.772,1.632,1.464,1.328,
     .  1.206,1.089,0.968,0.830,0.648,0.506,0.359,0.314,0.286,0.197,
     .  0.138,0.127/
      data mfe /0.007,0.015,0.018,0.024,0.046,0.103,0.149,0.222,0.298,
     .  0.380,0.471,0.576,0.703,0.872,0.991,1.088,1.112,1.093,1.028,
     .  1.014,0.990/
      data npa/7/
      data pamass /0.2026/
      data palogg /5.5941,6.2283,6.5230,6.6962,6.8549,6.9792,
     .  7.0418/

      Pi=acos(-1.e0)   !Pi
      G=6.674e-11 !N m^2 kg^-2  Gravitation constamt
      Msun=1.9891e30 !kg  mass of Sun
      Mearth=5.974e24 !kg mass of Earth
      Mjup=317.833e0*Mearth !kg  mass of Jupiter
      Rsun=696265.0e0*1000.0e0 !m  radius of Sun
      Rjup=142980.0e0*1000.0e0/2.0e0 !m  radius of Jupiter

      mmax=1.1 !maximum mass to plot and work with
 
      call pgopen('?')
      call pgask(.false.)
      call PGPAP ( 6.0 ,1.0) !square plot
      call pgslw(2) !thicker lines
      call pgsch(1.0) !bigger text

      call pgwindow(-6.0,log10(3.5),-2.2,log10(3.5))
            
      call pgbox('BCLNTS1',0.0,0,'BCLNTS1',0.0,0)
      call pglabel("Mass (M\d\(2281)\u)","Radius (R\d\(2281)\u)","")

cC     Carbon M-R relation
c      do 10 i=1,nc
c        rc(i)=log10(rc(i)/100.0)
c        mc(i)=log10(mc(i))
c 10   continue
c      call pgsci(2)
c      call pgline(nc,mc,rc)
c      call pgsci(1)

C     old Chandra relation  
c      do 12 i=1,nch
c        rch(i)=rch(i)/100.0
c 12   continue
c      call pgsci(4)
c      call pgline(nch,mch,rch)
c      call pgsci(1)
  
c      do 13 i=1,nfe
c        rfe(i)=log10(rfe(i)/100.0)
c        mfe(i)=log10(mfe(i))
c 13   continue
c      call pgsci(5)
c      call pgline(nfe,mfe,rfe)
c      call pgsci(1)
 
cC     New Chandra M-R relation
c      nch2=999
c      filename='chmassrad.dat'
c      open(unit=10,file=filename,status='old',err=901)       
c      do 15 i=1,nch2
c        read(10,*) mch2(i),rch2(i)
c        rch2(i)=log10(rch2(i))
c        mch2(i)=log10(mch2(i))
c 15   continue
c      call pgsci(6)
c      call pgline(nch2,mch2,rch2)
c      call pgsci(1)
c      close(10)
      
      do 25 j=1,5
        if(j.eq.1)filename='zp69.Mg.txt'
        if(j.eq.2)filename='zp69.He.txt'
        if(j.eq.3)filename='zp69.H.txt'
        if(j.eq.4)filename='zp69.C.txt'
        if(j.eq.5)filename='zp69.Fe.txt'
        open(unit=10,file=filename,status='old',err=901)
        read(10,*) dumc !first line is a comment
        i=1
 23     read(10,*,end=24) zpm(i),zpr(i)
            zpr(i)=log10(zpr(i)/100.0)
c        if(j.eq.3) then
c            write(6,*) zpm(i),zpr(i)
c            read(5,*)
c        endif
            i=i+1
        goto 23
 24     continue
        nzp=i-1
        call pgsci(8)
        call pgline(nzp,zpm,zpr)
        call pgsci(1)
 25   continue
      close(10)
      
      
      
C     Panei et al.MNRAS, 382, 779
      do 22 i=1,npa
        px(i)=log10(pamass)
        py(i)=log10(sqrt(pamass*1.989e33*6.67259e-8/
     .      10**palogg(i))/6.9599e10)
     
c        py(i)=log10(sqrt(6.9599e10*pamass*1.989e33/
c     .      (10**palogg(i)))/6.67259e-8)
c        write(6,*) px(i),py(i)
 22   continue
c      call pgsci(2)
c      call pgline(npa,px,py)
c      call pgpt(npa,px,py,21)
c      call pgsci(1)
        

cC     Simple R~M^-1/3 relation
c      dm=mmax/real(nmax)
c      do 11 i=1,nmax
c        msimp(i)=real(i)*dm
c        rsimp(i)=0.0084*msimp(i)**(-1.0/3.0)
c        msimp(i)=log10(msimp(i))
c        rsimp(i)=log10(rsimp(i))
c 11   continue
c      nsimp=nmax      
c      call pgsci(3)
c      call pgline(nsimp,msimp,rsimp)
c      call pgsci(1)
c      do 14 i=1,nmax
c        msimp(i)=real(i)*dm
c        rsimp(i)=0.01323*msimp(i)**(-1.0/3.0)
c        msimp(i)=log10(msimp(i))
c        rsimp(i)=log10(rsimp(i))
c 14   continue
c      nsimp=nmax      
c      call pgsci(3)
c      call pgline(nsimp,msimp,rsimp)
c      call pgsci(1)    

      call pgslw(1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
C     White Dwarf data
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Procyon B : Provencal et al. 1998
      mass(1)    =0.604 !Msun
      masserr(1) =0.018
      rad(1)     =0.0096 !Rsun
      raderr(1)  =0.0004
C     Sirius B : Provencal et al. 1998
      mass(2)    =1.000 !Msun
      masserr(2) =0.016
      rad(2)     =0.0084 !Rsun
      raderr(2)  =0.0002
C     40 Eri B
      mass(3)    =0.501 !Msun
      masserr(3) =0.011
      rad(3)     =0.0136 !Rsun
      raderr(3)  =0.0002
C     CD-38 10980
      mass(4)    =0.74 !Msun
      masserr(4) =0.04 
      rad(4)     =0.01245 !Rsun
      raderr(4)  =0.0004
C     W485A
      mass(5)    =0.59 !Msuin
      masserr(5) =0.04
      rad(5)     =0.0150 !Rsun
      raderr(5)  =0.001
C     L268-92
      mass(6)    =0.70
      masserr(6) =0.12
      rad(6)     =0.0149
      raderr(6)  =0.001
C     L481-60
      mass(7)    =0.53
      masserr(7) =0.05
      rad(7)     =0.0120
      raderr(7)  =0.0004
C     G154-B5B
      mass(8)    =0.46
      masserr(8) =0.08
      rad(8)     =0.013
      raderr(8)  =0.002
C     G181-B5B
      mass(9)    =0.50
      masserr(9) =0.05
      rad(9)     =0.011
      raderr(9)  =0.001
C     G156-64
      mass(10)   =0.59
      masserr(10)=0.06
      rad(10)    =0.0110
      raderr(10) =0.0010
C     NLTT 11748 Kwaka & Vennes 2009
      mass(11)   =0.167 !Msun
      masserr(11)=0.005
      rad(11)    =0.0540 !need to be checked
      raderr(11) =0.003 
C     PSR J1012+5307 van Kerkwuk 1996
      mass(12)   =0.16 !Msun
      masserr(12)=0.02
      rad(12)    =0.028 !logg=6.75+/-0.07 cgs (check)
      raderr(12) =0.005
C     PSR J1911-5958A Bassa 2006
      mass(13)   =0.18 !Msun
      masserr(13)=0.02
      rad(13)    =0.043 !logg=6.75+/-0.07 cgs
      raderr(13) =0.009
C     SDSS J123410.37-022802.9 Liebert 2004
      mass(14)   =0.185
      masserr(14)=0.005
      rad(14)    =0.046 !logg=6.38+/-0.05 (check)
      raderr(14) =0.009
      
      npt=14
      
      do 17 i=1,npt
        masserr(i)=log10(mass(i)+masserr(i))-log10(mass(i))
c        write(6,*) log10(mass(i)+masserr(i)),log10(mass(i))
        mass(i)=log10(mass(i))
        raderr(i)=log10(rad(i)+raderr(i))-log10(rad(i))
        rad(i)=log10(rad(i))
 17   continue

      call pgsch(1.0)
      call pgpt(npt,mass,rad,4)
      call pgerrb(5,npt,mass,rad,masserr,1.0)
      call pgerrb(6,npt,mass,rad,raderr,1.0)
      call pgsch(1.0)
      
C     Edmond, Gilliland object
      call pgpt1(log10(0.17),log10(0.105),4)
      
      filename="lowmassrad.dat"
      open(unit=10,file=filename,status='old',err=901)

 20   read(10,*,end=21) dumc,sma,smaerr,sra,sraerr
C       Stars
        smaerr=log10(sma+smaerr)-log10(sma)
        sma=log10(sma)
        sraerr=log10(sra+sraerr)-log10(sra)
        sra=log10(sra)
        call pgsci(6)
        call pgpt1(sma,sra,13)
        call pgerr1(5,sma,sra,smaerr,1.0)
        call pgerr1(6,sma,sra,sraerr,1.0)
      goto 20
 21   continue  
      close(10)
      call pgsci(1) !put colour back
      
      
      filename="exoplanets_massrad.csv"
      open(unit=10,file=filename,status='old',err=901)   
      read(10,500) dumc
 500  format(A1)
 
      i=0
 18   read(10,*,end=19) dumc,sma,smaerr,sra,sraerr,teff,pma,pmaerr,pra,
     . praerr
        i=i+1
C       Stars
        smaerr=log10(sma+smaerr)-log10(sma)
        sma=log10(sma)
        sraerr=log10(sra+sraerr)-log10(sra)
        sra=log10(sra)
        call pgsci(5)
        if(dumc(1:3).eq."Kep") call pgsci(4)
        call pgpt1(sma,sra,13)
        call pgerr1(5,sma,sra,smaerr,1.0)
        call pgerr1(6,sma,sra,sraerr,1.0)
C       Planets
        pma=pma*Mjup/Msun
        pmaerr=pmaerr*Mjup/Msun
        pra=pra*Rjup/Rsun
        praerr=praerr*Rjup/Rsun
        pmaerr=log10(pma+pmaerr)-log10(pma)
        pma=log10(pma)
        praerr=log10(pra+praerr)-log10(pra)
        pra=log10(pra)
        call pgsci(3)
        if(dumc(1:3).eq."Kep") call pgsci(4)
        call pgpt1(pma,pra,5)
        call pgerr1(5,pma,pra,pmaerr,1.0)
        call pgerr1(6,pma,pra,praerr,1.0)
c        write(6,*) pma,pra
      goto 18
 19   continue  
      close(10)
      call pgsci(1) !put colour back
      
      


CCCCCCCCCCCCCCCCCCCCCC      
C     KOI-74
      kmass(1)=   211*Mjup/Msun!18.6*Mjup/Msun!34.0*Mjup/Msun !Msun
      kmasserr(1)=10*Mjup/Msun!15.4*Mjup/Msun!2.7 *Mjup/Msun
      krad(1)=    0.383*Rjup/Rsun  !Rsun
      kraderr(1)= 0.010*Rjup/Rsun
C     Hoststar
      kmass(3)=   2.22 !Msun
      kmasserr(3)=0.12
      krad(3)=    1.899  !Rsun
      kraderr(3)= 0.05
      
C     KOI-81
      kmass(2)=   222.14*Mjup/Msun !guess
      kmasserr(2)=32.00*Mjup/Msun
      krad(2)=    1.119*Rjup/Rsun  !Rsun
      kraderr(2)= 0.0060*Rjup/Rsun
C     Hoststar
      kmass(4)=   2.71 !Msun
      kmasserr(4)=0.15
      krad(4)=    2.93  !Rsun
      kraderr(4)= 0.14
      
            
      do 16 i=1,4
        kmasserr1(i)=log10(kmass(i)+kmasserr(i))-log10(kmass(i))
        kmasserr2(i)=log10(kmass(i))-log10(kmass(i)-kmasserr(i))
        kmass(i)=log10(kmass(i))
        kraderr(i)=log10(krad(i)+kraderr(i))-log10(krad(i))
        krad(i)=log10(krad(i))
 16   continue
      
      call pgpt(4,kmass,krad,12)
      call pgerrb(1,4,kmass,krad,kmasserr1,1.0)
      call pgerrb(3,4,kmass,krad,kmasserr2,1.0)
      call pgerrb(6,4,kmass,krad,kraderr,1.0)
          
c      call pgsch(0.8)
c      call PGPTXT(-1.7, -1.35, 0.0, 0.0, 'KOI-74b')
c      call PGPTXT(-0.88, -0.88, 0.0, 0.0, 'KOI-81b')
c      call pgsch(1.0)
      
C     KOI-13 - dilution 0.09%
      kmass(1)=   0.052*Mjup/Msun!5.98*Mjup/Msun
      kmasserr(1)=0.033*Mjup/Msun!0.02*Mjup/Msun
      krad(1)=    0.63*Rjup/Rsun!1.21*Rjup/Rsun  !Rsun
      kraderr(1)= 0.04*Rjup/Rsun
C     KOI-13 - dilution 41%
      kmass(2)=   1.68*Mjup/Msun
      kmasserr(2)=0.1*Mjup/Msun
      krad(2)=    2.12*Rjup/Rsun
      kraderr(2)= 0.1*Rjup/Rsun

      do 26 i=1,2
        kmasserr1(i)=log10(kmass(i)+kmasserr(i))-log10(kmass(i))
        kmasserr2(i)=log10(kmass(i))-log10(kmass(i)-kmasserr(i))
        kmass(i)=log10(kmass(i))
        kraderr(i)=log10(krad(i)+kraderr(i))-log10(krad(i))
        krad(i)=log10(krad(i))
 26   continue
      
c      call pgpt(2,kmass,krad,12)
c      call pgerrb(1,2,kmass,krad,kmasserr1,1.0)
c      call pgerrb(3,2,kmass,krad,kmasserr2,1.0)
c      call pgerrb(6,2,kmass,krad,kraderr,1.0)
      
      call pgsci(2)
      call pgpt1(log10(Mearth/Msun),log10(0.009155),8)
      call pgpt1(log10(2.859e-4),log10(0.08241),11)
      call pgpt1(log10(Mjup/Msun),log10(Rjup/Rsun),11)
      call pgpt1(log10(5.1689e-5),log10(0.0353),11)
      call pgpt1(log10(4.3564e-5),log10(0.03633),11) !Uranus
      call pgpt1(0.0,0.0,9)
      
CCCCCCCCCCCCCCCCCCCCCCC      
      
      call pgclos()
      
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 999  end