      program koi701track
      implicit none
      integer nunit,nmax,i,dumi,npt,j,k,n1,nlast
      parameter(nmax=151)
      real Tsun, Pi, Rsun, Msun,dumr,logL,Teff(nmax),rhostar(nmax),
     .  trad,Ms,Rs,px(7),py(7),errxp(3),errxm(3),erryp(3),errym(3),
     .  M(7),logg(nmax),G,age(nmax),x,y,rbb(4),oldage,newage,age1,Teff1,
     .  logg1,Teff2,logg2,age2,page,dage,pTeff,plogg,aotu
      character*80 filename(7),dumc,text,tfile

      Tsun=5781.6e0 !Solar temperature (K)
      Pi=acos(-1.e0)!define Pi and 2*Pi
      Rsun=696265.0e0*1000.0e0 !m  radius of Sun
      Msun=1.9891e30 !kg  mass of Sun
      G=6.67259E-11
      aotu=14.0 !age of the Universe
      
      call pgopen('?')
      call pgpage()
      call PGPAP ( 6.0 ,1.0)  
      call pgvport(0.15,0.85,0.2,0.9)
c      call PGENV ( -0.5 , 0.5 , -0.5 , 0.5 , 1 , -2 ) 
      rbb(1)=12000.0
      rbb(2)=8000.0
      rbb(3)=4.5
      rbb(4)=3.0
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
c      call pglabel("T\deff\u (K)",
c     .  "\(0643) (g/cm\u3\d)","")

      call pglabel("T\deff\u (K)","log(g)","")
      
      dage=0.02 !delta age in Gyr
      call pgslw(1)
      call pgsch(1.0)
      call pgsci(15)
      do 16 i=200,0!349
        nunit=11
        write(tfile,501) "m",i,".track2"
 501    format(A1,I3,A7)
 
        ms=real(i)/100.0*Msun
        open(unit=nunit,file=tfile,status='old',err=902)
        
        do 19 j=1,11
            read(nunit,*) dumc
 19     continue
        
        
        read(nunit,*,end=17) dumi,age1,Teff1,logL
        Teff1=10.0**(Teff1)
        trad=sqrt( (Tsun/Teff1)**4.0e0*10.0e0**logL )
        Rs=trad*Rsun
        logg1=log10(G*Ms/(Rs*Rs))+2.0
        page=age1
 18     read(nunit,*,end=17) dumi,age2,Teff2,logL
c        if(age2.gt.aotu) goto 17 !don't plot tracks beyond 14 Gyr)
            Teff2=10.0**(Teff2)
            trad=sqrt( (Tsun/Teff2)**4.0e0*10.0e0**logL )
            Rs=trad*Rsun
            logg2=log10(G*Ms/(Rs*Rs))+2.0
c            write(0,*) Teff2,logg2
c            page=age1
            k=0
            do 20 while(page.le.age2)
C               do a linear interpolation
                pTeff=Teff1+(Teff1-Teff2)*(page-age1)/(age2-age1)  
                plogg=logg1+(logg1-logg2)*(page-age1)/(age2-age1)         
                call pgpt1(pTeff,plogg,1)
                k=k+1
                page=page+dage !increase age
 20         continue
c            page=page-dage
            if(k.gt.0)then
                Teff1=Teff2
                logg1=logg2
                age1=age2
            endif
        goto 18
 17     continue       
        close(nunit)
c        read(5,*)
 16   continue
      call pgslw(3)
      call pgsch(1.5)
      call pgsci(1)
      
      
      errxp(1)=250.0
      errxm(1)=250.0
      errxp(2)=500.0
      errxm(2)=500.0
      errxp(3)=750.0
      errxm(3)=750.0
      erryp(1)=0.0092
      errym(1)=0.0092
      erryp(2)=0.0184
      errym(2)=0.0184
      erryp(3)=0.0276
      errym(3)=0.0276
      
      do 15 j=1,3
        px(1)=9974.0-errxm(j)
        px(2)=9974.0+errxp(j)
        px(3)=px(2)
        px(4)=px(1)
        px(5)=px(1)

        py(1)=4.0608-errym(j)
        py(2)=py(1)
        py(3)=4.0608+erryp(j)
        py(4)=py(3)
        py(5)=py(1)
      
c      do 14 i=1,5
c        write(0,*) px(i),py(i)     
c 14   continue
      
        call pgline(5,px,py)

 15   continue
     
      nunit=10
      filename(1)="m200.track2"
      filename(2)="m230.track2"
      filename(3)="m260.track2"
      filename(4)="m290.track2"
      filename(5)="m320.track2"
      filename(6)="m350.track2"
      filename(7)="m380.track2"
      M(1)=2.0
      M(2)=2.3
      M(3)=2.6
      M(4)=2.9
      M(5)=3.2
      M(6)=3.5
      M(7)=3.8
      
      do 13 j=1,7
        open(unit=nunit,file=filename(j),status='old',err=901)
        call pgsci(j+1)
      
        do 10 i=1,11
            read(nunit,*) dumc
 10     continue
      
        n1=0
        oldage=0.0
        nlast=0
       
        Ms=M(j)*Msun
        i=1
 11     read(nunit,*,end=12) dumi,age(i),Teff(i),logL
            Teff(i)=10.0**(Teff(i))
            trad=sqrt( (Tsun/Teff(i))**4.0e0*10.0e0**logL )
            Rs=trad*Rsun
            rhostar(i)=Ms/(4.0d0/3.0d0*pi*Rs*Rs*Rs)/1000.0d0
c            write(0,*) i,Teff(i),rhostar(i)
            logg(i)=log10(G*Ms/(Rs*Rs))+2.0
c            write(0,*) Teff(i),logg(i),M(j),trad
            call pgsch(1.0)
            write(text,500) floor(age(i)*10.0)/10.0
 500        format(F4.1)
            x=teff(i)
            y=logg(i)
            newage=floor(age(i)*10.0)/10.0
            if((x.lt.rbb(1)).and.(x.gt.rbb(2)+60).and.
     .        (y.lt.rbb(3)).and.(y.gt.rbb(4)))then
                if((n1.eq.0).or.(newage.gt.oldage+0.1))then
                    if((n1.eq.0).or.(nlast.ge.9))then
                        call pgtext(x,y,text)
                        call pgpt1(x,y,2)
                        n1=1
                        oldage=newage
                        nlast=0
                    endif
                endif
            endif
            nlast=nlast+1           

            call pgsch(1.5)

c            if(mod(i+j,8).eq.0)then
c                write(text,502) age(i)
c 502            format(F4.1)
c                x=teff(i)
c                y=logg(i)
c                call pgsch(1.0)
c                if((x.lt.rbb(1)).and.(x.gt.rbb(2)+50).and.
c     .              (y.lt.rbb(3)).and.(y.gt.rbb(4)))then
c                    call pgtext(x,y,text)
c                    call pgpt1(x,y,2)
c                endif
c                call pgsch(1.5)
c            endif

            i=i+1
        goto 11     
 12     continue
        npt=i-1
      
        close(nunit)
      
        call pgline(npt,Teff,logg)
 
 13   continue    
      call pgsci(1)
      

      call pgclos()

      goto 999
 901  write(0,*) "Cannot open ",filename(j)
      goto 999
 902  write(0,*) "Cannot open ",tfile
      goto 999
 999  end
