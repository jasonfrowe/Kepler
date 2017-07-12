      program TemperaturevsAlbedo
      implicit none
      integer nmax,i,j,nstar,k
      parameter(nmax=10000,nstar=4)
      real Pi,G,Per(nstar),Rad(nstar),Ms(nstar),Mp(nstar),Ts(nstar),
     .  ab(nstar),aberr(nstar),Msun,Mearth,Mjup,Rsun,lx(200),ly(200),
     .  Psec,M1,M2,Rs,asemi,Rsd2a
      character*5 kname(4)

      call pgopen('?')
      call PGPAP ( 4.0 ,1.0)  
c      call PGENV ( -0.5 , 0.5 , -0.5 , 0.5 , 1 , -2 )  
      call pgwindow(0.0,0.50,1400.0,2400.0)
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Bond Albedo","T\deq\u","")


      Pi=acos(-1.0e0)   !Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constamt
      
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Mjup=317.833*Mearth !kg  mass of Jupiter
      Rsun=696265.0*1000.0 !m  radius of Sun      

C     KOI-1
      Per(1)=2.4706131924
      Rad(1)=1.4719
      Ms(1)=0.99945
      Mp(1)=1.5
      Ts(1)=5250.0
      ab(1)=0.04
      aberr(1)=0.02
      kname="KOI-1"


C     KOI-18      
c      Per(1)=3.548461
c      Rad(1)=2.066
c      Ms(1)=1.448
c      Mp(1)=2.164
c      Ts(1)=6297.0
c      ab(1)=0.093
c      aberr(1)=0.041
c      kname(1)="koi18"
C     KOI-97
      Per(2)=4.885525
      Rad(2)=1.855
      Ms(2)=1.371
      Mp(2)=0.426
      Ts(2)=6000.0
      ab(2)=0.142
      aberr(2)=0.059
      kname(2)="koi97"
C     KOI-7
      Per(3)=3.21339
      Rad(3)=1.588
      Ms(3)=1.250
      Mp(3)=0.076
      Ts(3)=5857.0
      ab(3)=0.152
      aberr(3)=0.437
      kname(3)="koi07"
C     KOI-17
      Per(4)=3.234719
      Rad(4)=1.386
      Ms(4)=1.217
      Mp(4)=0.649
      Ts(4)=5826.0
      ab(4)=0.034
      aberr(4)=0.054
      kname(4)="koi17"   
            
      j=1
      do 10 k=1,1
        ab(k)=ab(k)*3.0/2.0
        aberr(k)=aberr(k)*3.0/2.0
        Psec=Per(k)*24.0*60.0*60.0 !sec period of planet in sec
        M1=Ms(k)*Msun!1.083*Msun !kg  mass of star
        M2=Mp(k)*Mjup!0.69*Mjup !kg   mass of planet
        Rs=Rad(k)*Rsun!1.083*Rsun
        asemi=(Psec*Psec*G*(M1+M2)/(4.0*Pi*Pi))**(1.0/3.0)
      
        Rsd2a=sqrt(Rs/(2.0*asemi))
      
c        lx(1)=0.0
c        ly(1)=Ts(k)*Rsd2a*(1.0*(1- (0.0) ))**0.25
c        lx(2)=ab(k)+aberr(k)
c        ly(2)=Ts(k)*Rsd2a*(1.0*(1- (ab(k)+aberr(k)) ))**0.25
c        lx(3)=ab(k)+aberr(k)
c        ly(3)=Ts(k)*Rsd2a*(2.0*(1- (ab(k)+aberr(k)) ))**0.25
c        lx(4)=0.0
c        ly(4)=Ts(k)*Rsd2a*(2.0*(1- (0.0) ))**0.25

        do 11 i=1,100
            lx(i)=(ab(k)+aberr(k))*real(i-1)/100.0
            ly(i)=Ts(k)*Rsd2a*(1.0*(1- (lx(i)) ))**0.25
            lx(i+100)=(ab(k)+aberr(k))*real(100-i)/100.0
            ly(i+100)=Ts(k)*Rsd2a*(2.0*(1- (lx(i+100)) ))**0.25
c            write(6,*) i,lx(i),ly(i)
c            write(6,*) i+100,lx(i+100),ly(i+100)
c            read(5,*)
 11     continue
        
c        write(6,500) kname(k),ly
c 500    format(A5,1X,4(F8.3,1X))
        write(6,*) kname(k), Ts(k)*Rsd2a*(1.0*(1- (0.1) ))**0.25
        
        j=j+1
        call pgsci(j)
        call pgslw(1)
        call pgsls(1)
        call pgsfs(j-1)
        call pgpoly(200,lx,ly)
        call pgsfs(1)
c        lx(1)=ab(k)
c        ly(1)=Ts(k)*Rsd2a*(1.0*(1- (ab(k)) ))**0.25
c        lx(2)=ab(k)
c        ly(2)=Ts(k)*Rsd2a*(2.0*(1- (ab(k)) ))**0.25
c        call pgline(2,lx,ly)
 10   continue
      
      call pgclos()
      
      end