      program koi20track
      implicit none
      integer nunit,nmax,i,dumi,npt,j
      parameter(nmax=151)
      real Tsun, Pi, Rsun, Msun,dumr,logL,Teff(nmax),rhostar(nmax),
     .  trad,Ms,Rs,px(5),py(5),errxp(3),errxm(3),erryp(3),errym(3),
     .  logg(nmax),G
      character*80 filename(5),dumc

      Tsun=5781.6e0 !Solar temperature (K)
      Pi=acos(-1.e0)!define Pi and 2*Pi
      Rsun=696265.0e0*1000.0e0 !m  radius of Sun
      Msun=1.9891e30 !kg  mass of Sun
      G=6.67259E-11
      
      call pgopen('?')
      call pgpage()
      call PGPAP ( 6.0 ,1.0)  
      call pgvport(0.15,0.85,0.2,0.9)
c      call PGENV ( -0.5 , 0.5 , -0.5 , 0.5 , 1 , -2 )  
      call pgwindow(7500.0,5000.0,5.0,2.0)
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("T\deff\u (K)",
     .  "\(0643) (g/cm\u3\d)","")
     
      nunit=10
      filename(1)="m13x65z04.track2"
      filename(2)="m16x65z04.track2"
      filename(3)="m19x65z04.track2"
      filename(4)="m22x65z04.track2"
      filename(5)="m25x65z04.track2"
      
      do 13 j=1,5
        open(unit=nunit,file=filename(j),status='old',err=901)
      
        do 10 i=1,11
            read(nunit,*) dumc
 10     continue
      
        Ms=1.1*Msun
        i=1
 11     read(nunit,*,end=12) dumi,dumr,Teff(i),logL
            Teff(i)=10.0**(Teff(i))
            trad=sqrt( (Tsun/Teff(i))**4.0e0*10.0e0**logL )
            Rs=trad*Rsun
            rhostar(i)=Ms/(4.0d0/3.0d0*pi*Rs*Rs*Rs)/1000.0d0
c            write(0,*) i,Teff(i),rhostar(i)
            logg(i)=log10(G*Ms/(Rs*Rs))+2.0
            i=i+1
        goto 11     
 12     continue
        npt=i-1
      
        close(nunit)
        call pgsci(j+1)
        call pgline(npt,Teff,logg)
 
 13   continue    
      call pgsci(1)
      
      errxp(1)=250.0
      errxm(1)=250.0
      errxp(2)=500.0
      errxm(2)=500.0
      errxp(3)=750.0
      errxm(3)=750.0
      erryp(1)=0.2
      errym(1)=0.2
      erryp(2)=0.4
      errym(2)=0.4
      erryp(3)=0.6
      errym(3)=0.6
      
      do 15 j=1,3
        px(1)=6300.0-errxm(j)
        px(2)=6300.0+errxp(j)
        px(3)=px(2)
        px(4)=px(1)
        px(5)=px(1)

        py(1)=3.6-errym(j)
        py(2)=py(1)
        py(3)=3.6+erryp(j)
        py(4)=py(3)
        py(5)=py(1)
      
c      do 14 i=1,5
c        write(0,*) px(i),py(i)     
c 14   continue
      
        call pgline(5,px,py)

15    continue

      call pgclos()

      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 999  end