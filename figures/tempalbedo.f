      program TemperaturevsAlbedo
      implicit none
      integer nmax,i,j
      parameter(nmax=10000)
      real ab(nmax),temp(nmax),Ts,Rs,asemi,f,Per,Psec,Pi,G,M1,M2,Msun,
     .   Mearth,Mjup,rnmax,Rsd2a,Rsun,lx(2),ly(2)
      character*80 label


      Pi=acos(-1.0e0)   !Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constamt
      Per=2.4706131924
      Psec=Per*24.0*60.0*60.0 !sec period of planet in sec
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Mjup=317.833*Mearth !kg  mass of Jupiter
      M1=0.99945*Msun!1.083*Msun !kg  mass of star
      M2=1.4719*Mjup!0.69*Mjup !kg   mass of planet
      asemi=(Psec*Psec*G*(M1+M2)/(4.0*Pi*Pi))**(1.0/3.0)
      Rsun=696265.0*1000.0 !m  radius of Sun
      Rs=0.96358*Rsun!1.083*Rsun
      Ts=5250.0!6000.0
      
      Rsd2a=sqrt(Rs/(2.0*asemi))
      rnmax=real(nmax)

      lx(1)=0.020
      lx(2)=0.024

      call pgopen('?')
      call PGPAP ( 6.0 ,1.0)  
      call PGENV ( -0.5 , 0.5 , -0.5 , 0.5 , 1 , -2 )  
      call pgwindow(0.0,0.03,1200.0,1700.0)
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTSV1',0.0,0)

      f=1.0
      ly(1)=1615.0
      ly(2)=1615.0
      do 11 j=1,3
         call pgsls(j)
         do 10 i=1,nmax
            ab(i)=real(i-1)/rnmax
            temp(i)=Ts*Rsd2a*(f*(1-ab(i)))**0.25
 10      continue
         call pgline(nmax,ab,temp)
         call pgline(2,lx,ly)
         write(label,500) "f = ",f
 500     format(A4,F3.1)
         call pgtext(0.025,ly(1),label)
         f=f+0.5
         ly(1)=ly(1)+25.0
         ly(2)=ly(2)+25.0
 11   continue    
      
      call pglabel("Bond Albedo","T\deq\u","")
      call pgslw(1)
c      lx(1)=0.0
c      lx(2)=1.0
c      ly(1)=1130.0
c      ly(2)=1130.0
c      call pgsls(2)
c      call pgline(2,lx,ly)
      lx(1)=0.009
      lx(2)=0.009
      ly(1)=1200.0
      ly(2)=1700.0
      call pgline(2,lx,ly)
      call pgsls(1)
      call pgsfs(3)
      call pgrect(0.0,lx(2),ly(1),ly(2))
      call pgsfs(1)
      call pgclos()
      
      end