      program nfitplots
      implicit none
      integer nmax,iargc,nunit,i,npt,nbinmax,j,now(3),seed,nabin,
     .  nfp,kid,fpflag,nKepler,k
      parameter(nmax=10000,nbinmax=500,nabin=10,nKepler=491)
      integer flag(nmax),fp(nmax)
      real mstar(nmax),merr(nmax),rstar(nmax),rerr(nmax),Teff(nmax),
     .  logg(nmax),FeH(nmax),Per(nmax),rhofit(nmax),rhostar(nmax),
     .  rdr(nmax),rdrerr(nmax),rp(nmax),rperr(nmax),bn(nmax),ba(nmax),
     .  occ(nmax),occerr(nmax),rbb(4),px(nmax),py(nmax),x,y,
     .  bdatax(nbinmax),bdatay(nbinmax),asemi(nmax),Teq(nmax),td(nmax),
     .  tdur1,tdur2,transitdur,dumr,ran2,ranb,ave,adev,sdev,var,skew,
     .  curt,mean,radearth(nmax),ab(nmax),Tbin(nabin),Abin(nabin),
     .  esum(nabin),Asum(nabin),Abinerr(nabin),Aerr,koi(nmax),Pi,
     .  koinum(nmax),aberr(nmax),Fin(nmax),Fin1,Fin2,AU,Rsun,sigma,Fsun,
     .  Keplerlam(nKepler),Keplerpass(nKepler),st_temp,st_flux,planck,
     .  lx(nmax),ly(nmax),Rjup,dist,Dp,Abcor(nmax)
      character*80 filename,koichar,title,fplist
      
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      sigma=5.67d-8 !stephan-boltzman const.
      AU=1.49598d11 !astronomical unit
      Pi=acos(-1.d0)   !Pi
     
      nunit=10
      fplist="koi_FPs.20120229.csv"
      i=1
      open(unit=nunit,file=fplist,status='old',err=903)
 5    read(nunit,*,end=6) koi(i),kid,dumr,fp(i)
        if(fp(i).eq.1)then
c            write(0,*) koi(i),kid,fp(i)
            i=i+1
        endif
      goto 5
 6    continue
      nfp=i-1
      write(0,*) "nFP: ",nfp
     
      call itime(now)
      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      dumr=ran2(-seed)
           
      if(iargc().lt.1) goto 901
      
      call getarg(1,filename)
      open(unit=nunit,file=filename,status='old',err=901)
      
      k=1
      i=1
 10   read(nunit,*,err=904,end=11) koichar,mstar(i),merr(i),rstar(i),
     .  rerr(i),
     .  Teff(i),logg(i),FeH(i),Per(i),rhofit(i),rhostar(i),rdr(i),
     .  rdrerr(i),rp(i),rperr(i),bn(i),ba(i),occ(i),occerr(i),flag(i)
c        write(0,*) koichar,mstar(i)
        
        fpflag=0
        read(koichar(2:9),*) koinum(i)
c       write(0,*) koinum-int(koinum)*100)
c        if(int( (koinum-int(koinum))*100 ).ne.2) fpflag=1
        do 25 j=1,nfp
            if(int(koinum(i)).eq.int(koi(j)))then
                fpflag=1
            endif
 25     continue
cc        if(occerr(i).le.0)fpflag=1
cc        if(occ(i)/occerr(i).lt.3.0) fpflag=1
c        fpflag=0
c        bn(i)=bn(i)*bn(i)
        if(fpflag.eq.0) then
c            write(6,*) koichar(1:9)
            i=i+1
        endif
        k=k+1
        goto 10
 11   continue
      npt=i-1
      write(0,*) "Points read ",npt
      close(nunit)
      
      call pgopen('?')
c      call pgpage()
      call PGPAP ( 6.0 ,1.0)  
      call pgvport(0.15,0.85,0.2,0.9)
      call pgpage()
      
      rbb(1)=log10(0.5)
      rbb(2)=log10(600.0)
      rbb(3)=-4.0
      rbb(4)= 4.0
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCLNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Period (d)","Stellar Density (g/cm\u3\d)","")
      
      j=0
      call pgbbuf()
      mean=0.0
      do 12 i=1,npt
        x=log10(Per(i))
        y=rhostar(i)-rhofit(i)
        call pgsci(flag(i)+1)
        if((x.ge.rbb(1)).and.(x.le.rbb(2)).and.(y.ge.rbb(3)).and.
     .      (y.le.rbb(4))) then
            call pgpt1(x,y,17)
        endif
        if(abs(y).lt.4.0)then
            j=j+1
            px(j)=y
            mean=mean+y
        endif
 12   continue
      call pgebuf()
      call pgsci(1)
      mean=mean/real(j)
      write(0,*) "Mean drhostar :",mean
      
      title='delta rhostar'
      call histogram(30,j,px,nbinmax,bdatax,bdatay,title)
      
      title='impact parameter'
      j=0
      do 17 i=1,npt
        if(bn(i)*bn(i).le.1.0)then
            j=j+1
            px(j)=bn(i)*bn(i)
        endif
 17   continue
      call histogram(30,j,px,nbinmax,bdatax,bdatay,title)
      
      title='impact parameter'
      j=0
      do 18 i=1,npt
        if(ba(i)*ba(i).le.1.0)then
            j=j+1
            px(j)=ba(i)*ba(i)
        endif
 18   continue     
      call histogram(30,j,px,nbinmax,bdatax,bdatay,title)
      
      j=0
      do 16 i=1,npt
        x=ranb(Per(i),rhofit(i),seed)
        if((x.ge.0.0).and.(x.le.1.0))then
            j=j+1
            px(j)=x
        endif
 16   continue
      call histogram(30,j,px,nbinmax,bdatax,bdatay,title)
      
      
      do 13 i=1,npt
        asemi(i)=(per(i)*per(i)*10.0**logg(i)*rstar(i)*rstar(i)*
     .      9.17e23)**(1.0/3.0)/1.49598e11
        Teq(i)=Teff(i)*0.0682*(rstar(i)/(2.0*asemi(i)))**(1.0/2.0)*
     .      0.9146912
        radearth(i)=rp(i)*6378.1*1000.0
        ab(i)=occ(i)*(1.49598e11*asemi(i)/radearth(i))**2.0*1.0e-6
        aberr(i)=(occ(i)+occerr(i))*(1.49598e11*asemi(i)/
     .      radearth(i))**2.0*1.0e-6
        aberr(i)=aberr(i)-ab(i)
 13   continue
 
      do 22 i=1,nabin
        esum(i)=0.0
        Tbin(i)=0.0
        Abin(i)=0.0
        Asum(i)=0.0
        Tbin(i)=500.0+(3000.0-500.0)/real(nabin)*real(i-1)+125.0
 22   continue
      
      do 21 i=1,npt
c        j=int(Teq(i)/100.0-4.0+0.5)
        j=int((Teq(i)-500.0)/((3000.0-500.0)/real(nabin))+1.0)
        if((j.ge.1).and.(j.le.nabin).and.(ab(i).lt.1.0).and.
     .    (ab(i).gt.-1.0).and.(rp(i).gt.7.0).and.(occerr(i).gt.0.0)
     .    .and.(rp(i).lt.15.0))then
            Aerr=occerr(i)*(1.49598e11*asemi(i)/radearth(i))**2.0*1.0e-6
            esum(j)=esum(j)+1.0/Aerr
            Abin(j)=Abin(j)+ab(i)/Aerr
            Asum(j)=Asum(j)+1.0
        endif
 21   continue
 
      do 23 i=1,nabin
        if(esum(i).gt.0.0)then
            Abin(i)=Abin(i)/esum(i)
            Abinerr(i)=sqrt(Asum(i))/esum(i)
        endif
c        write(0,*) Tbin(i),Abin(i),Abinerr(i)
 23   continue
      
      call pgpage()
      rbb(1)=500.0
      rbb(2)=2000.0
      rbb(3)=-0.0
      rbb(4)= 200.0
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Teq (K)","Occult (ppm)","")
      
      call pgbbuf()
      do 14 i=1,npt
        if((teq(i).gt.rbb(1)).and.(teq(i).lt.rbb(2)))then
            if(rp(i).gt.7.0)then
                call pgsci(2)
            elseif(rp(i).gt.2.0)then
                call pgsci(3)
            else
                call pgsci(1)
            endif
            if(occ(i)/occerr(i).gt.2.0)then
c                write(6,*) koinum(i),occ(i),occerr(i)
                call pgpt1(teq(i),occ(i),17)
                call pgerr1(6,teq(i),occ(i),occerr(i),0.0)
            endif
        endif
 14   continue     
      call pgebuf()
      call pgsci(1)

      call pgpage()
      rbb(1)=500.0
      rbb(2)=3000.0
      rbb(3)=-0.1
      rbb(4)= 1.0
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Teq (K) (Ab=0.1)","Ag","")

c      call pgsci(3)
c      do 24 i=1,nabin
c        if(Abin(i).ne.0.0)then
c            call pgpt1(Tbin(i),Abin(i),17)
c            call PGERR1(1,tbin(i),abin(i),abinerr(i),1.0)
c        endif
c 24   continue
c      call pgsci(1)

      
      call pgbbuf()
      do 20 i=1,npt
        if((teq(i).gt.rbb(1)).and.(teq(i).lt.rbb(2)).and.
     .     (ab(i).gt.rbb(3)).and.(ab(i).lt.rbb(4)))then
            if(rp(i).gt.7.0)then
                call pgsci(2)
            elseif(rp(i).gt.2.0)then
                call pgsci(3)
            else
                call pgsci(1)
            endif
c            if((occ(i)/occerr(i).gt.2.0).and.(ab(i).lt.1.0))then
            if((aberr(i).lt.0.3).and.(ab(i).lt.1.0).and.
     .         (aberr(i).ne.0))then
                call pgpt1(teq(i),ab(i),17)
                call pgerr1(6,teq(i),ab(i),aberr(i),0.0)
c                write(6,500) koinum(i),int(teq(i)),ab(i),aberr(i)
 500            format(F7.2,1X,I4,1X,F5.2,1X,F5.2)
            endif
        endif
 20   continue     
      call pgebuf()
      call pgsci(1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      call pgpage()
      rbb(1)=1.0
      rbb(2)=5.0
      rbb(3)=-0.1
      rbb(4)= 1.0
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BLNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("F*/Fearth","Ag","")

      do 26 i=1,npt
        Fin(i)=Rstar(i)**2.0*(Teff(i)/5770.0d0)**4.0*asemi(i)**-2.0
        if((aberr(i).lt.0.3).and.(ab(i).lt.1.0).and.
     .         (aberr(i).ne.0).and.(ab(i)+aberr(i).gt.0.0))then

            if(rp(i).gt.7.0)then
                call pgsci(2)
            elseif(rp(i).gt.2.0)then
                call pgsci(3)
            else
                call pgsci(1)
            endif

            write(6,*) koinum(i),Teq(i),Ab(i)
            call pgpt1(log10(Fin(i)),Ab(i),17)
            call pgerr1(6,log10(Fin(i)),Ab(i),aberr(i),0.0)
        endif
 26   continue
      call pgsci(1)
 
      dp=AU*sqrt(1/10.0)
      Fin1=5770.0*sqrt(Rsun/(2.0*dp))*(1.0*(1-0.3)) !get Teq at bounds 
      dp=AU*sqrt(1/100000.0)
      Fin2=5770.0*sqrt(Rsun/(2.0*dp))*(1.0*(1-0.3))     

c      write(6,*) "F:",Fin1,Fin2
      call pgwindow(real(Fin1),real(Fin2),rbb(3),rbb(4))
      CALL PGBOX('CMTS1',0.0,0,'',0.0,0)
      
      call readKeplerbandpass(nKepler,Keplerlam,Keplerpass)
      st_temp=5777.0
      st_flux=Planck(st_temp,nKepler,Keplerlam,Keplerpass)
      st_flux=st_flux*4.0*Pi*6.9599d10*6.9599d10
      
c      Rjup=142980.0*1000.0*100.0/2.0 !cm  radius of Jupiter
      Rjup=69173.0*1000.0*100.0
           
      do 27 i=1,nmax
        lx(i)=Fin1+(Fin2-Fin1)/real(nmax)*real(i-1)
        ly(i)=Planck(lx(i),nKepler,Keplerlam,Keplerpass)
        ly(i)=ly(i)*4.0*Pi*(1.0*Rjup)**2.0/st_flux
        dist=Rsun/2.0*(5770.0/lx(i))**2.0*sqrt(1.0*(1.0-0.1))
        dist=dist*100.0d0 !change to cm
        ly(i)=ly(i)*(dist/(1.0*Rjup))**2.0
c        write(6,*) lx(i),ly(i),dist/(AU*100.0)
c        read(5,*)
 27   continue
c      call pgline(nmax,lx,ly)

      do 31 i=1,nmax
        lx(i)=Fin1+(Fin2-Fin1)/real(nmax)*real(i-1)
        ly(i)=Planck(lx(i),nKepler,Keplerlam,Keplerpass)
        ly(i)=ly(i)*4.0*Pi*(1.0*Rjup)**2.0/st_flux
        dist=Rsun/2.0*(5770.0/lx(i))**2.0*sqrt(2.0*(1.0-0.1))
        dist=dist*100.0d0 !change to cm
        ly(i)=ly(i)*(dist/(1.0*Rjup))**2.0
c        write(6,*) lx(i),ly(i),dist/(AU*100.0)
c        read(5,*)
 31   continue
c      call pgline(nmax,lx,ly)

      st_temp=4000.0
      st_flux=Planck(st_temp,nKepler,Keplerlam,Keplerpass)
      st_flux=st_flux*4.0*Pi*(0.5*6.9599d10)**2.0      
      do 28 i=1,nmax
        lx(i)=Fin1+(Fin2-Fin1)/real(nmax)*real(i-1)
        ly(i)=Planck(lx(i),nKepler,Keplerlam,Keplerpass)
        ly(i)=ly(i)*4.0*Pi*(Rjup/10.863)**2.0/st_flux
        dist=Rsun/2.0*(5770.0/lx(i))**2.0*sqrt(1.0*(1.0-0.1))
        dist=dist*100.0d0 !change to cm
        ly(i)=ly(i)*(dist/(Rjup/10.863))**2.0
c        write(6,*) lx(i),ly(i),dist/(AU*100.0)
c        read(5,*)
 28   continue
c      call pgline(nmax,lx,ly)
      
      st_temp=4000.0
      st_flux=Planck(st_temp,nKepler,Keplerlam,Keplerpass)
      st_flux=st_flux*4.0*Pi*(0.5*6.9599d10)**2.0      
      do 30 i=1,nmax
        lx(i)=Fin1+(Fin2-Fin1)/real(nmax)*real(i-1)
        ly(i)=Planck(lx(i),nKepler,Keplerlam,Keplerpass)
        ly(i)=ly(i)*4.0*Pi*(Rjup/10.863)**2.0/st_flux
        dist=Rsun/2.0*(5770.0/lx(i))**2.0*sqrt(2.0*(1.0-0.1))
        dist=dist*100.0d0 !change to cm
        ly(i)=ly(i)*(dist/(Rjup/10.863))**2.0
c        write(6,*) lx(i),ly(i),dist/(AU*100.0)
c        read(5,*)
 30   continue
c      call pgline(nmax,lx,ly)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Correct for thermal radiation
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      call pgpage()
      rbb(1)=1.0
      rbb(2)=5.0
      rbb(3)=-0.1
      rbb(4)= 1.0
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BLNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("F*/Fearth","Ag","")

      do 29 i=1,npt
c        Fin(i)=Rstar(i)**2.0*(Teff(i)/5770.0d0)**4.0*asemi(i)**-2.0

        st_temp=Teff(i)
        st_flux=Planck(st_temp,nKepler,Keplerlam,Keplerpass)
        st_flux=st_flux*4.0*Pi*(Rstar(i)*6.9599d10)**2.0 
        
        Abcor(i)=Planck(Teq(i),nKepler,Keplerlam,Keplerpass)
        Abcor(i)=Abcor(i)*4.0*Pi*(rp(i)/10.863*Rjup)**2.0/st_flux
        dist=Rsun/2.0*(5770.0/Teq(i))**2.0*sqrt(2.0*(1.0-0.1))
        dist=dist*100.0d0 !change to cm
        Abcor(i)=Abcor(i)*(dist/(rp(i)/10.863*Rjup))**2.0

        if((aberr(i).lt.0.3).and.(ab(i).lt.1.0).and.
     .         (aberr(i).ne.0).and.(ab(i)+aberr(i).gt.0.0))then

            if(rp(i).gt.7.0)then
                call pgsci(2)
            elseif(rp(i).gt.2.0)then
                call pgsci(3)
            else
                call pgsci(1)
            endif

            write(6,*) koinum(i),Ab(i),Abcor(i)
            call pgpt1(log10(Fin(i)),Ab(i)-Abcor(i),17)
            call pgerr1(6,log10(Fin(i)),Ab(i)-Abcor(i),aberr(i),0.0)
        endif
 29   continue
      call pgsci(1)

c      write(6,*) "F:",Fin1,Fin2
      call pgwindow(real(Fin1),real(Fin2),rbb(3),rbb(4))
      CALL PGBOX('CMTS1',0.0,0,'',0.0,0)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call pgpage()
      rbb(1)=200.0
      rbb(2)=2000.0
      rbb(3)=0.0
      rbb(4)=10.0
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Teq (K)","Radius (Rearth)","")
      
      call pgbbuf()
      do 19 i=1,npt
        if((teq(i).gt.rbb(1)).and.(teq(i).lt.rbb(2)).and.
     .     (rp(i).gt.rbb(3)).and.(rp(i).lt.rbb(4)))then
            if(rp(i).gt.7.0)then
                call pgsci(2)
            elseif(rp(i).gt.2.0)then
                call pgsci(3)
            else
                call pgsci(1)
            endif
            call pgpt1(teq(i),rp(i),17)
        endif
 19   continue     
      call pgebuf()      
      call pgsci(1)

      call pgpage()
      rbb(1)=log10(0.5)
      rbb(2)=log10(600.0)
      rbb(3)=0.0
      rbb(4)=2.0
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCLNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Period (d)","Td","")
      
      call pgbbuf()
      j=0
      do 15 i=1,npt
        tdur1=transitdur(Per(i),rhofit(i),rdr(i),bn(i))
        tdur2=transitdur(Per(i),rhostar(i)-0.4,rdr(i),0.0)!bn(i))
c        write(0,*) tdur1,tdur2
        if((tdur1/tdur2.gt.rbb(3)).and.(tdur1/tdur2.lt.rbb(4)))then
c            if(bn(i).gt.0.1)then
                x=log10(Per(i))
                y=tdur1/tdur2
                j=j+1
                td(j)=y
                if((x.gt.rbb(1)).and.(x.lt.rbb(2)).and.(y.gt.rbb(3))
     .              .and.(y.lt.rbb(4)))then
                    call pgpt1(x,y,17)
                endif
c            endif
        endif
 15   continue
      call pgebuf()
      write(0,*) "j: ",j
      
      title='Td'
      call histogram(30,j,td,nbinmax,bdatax,bdatay,title)
      call moment(td,j,ave,adev,sdev,var,skew,curt)
      write(0,*) "ave :",ave
      write(0,*) "adev:",adev
      write(0,*) "sdev:",sdev
      write(0,*) "var :",var
      write(0,*) "skew:",skew
      write(0,*) "curt:",curt
      
      call pgclos()
      goto 999
 901  write(0,*) "Usage: nfitplots nfitpars"
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 903  write(0,*) "Cannot open ",fplist
      goto 999
 904  write(0,*) "Error on line:",k
      goto 999
 999  end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function planck(T,npass,lam,pass)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     I believe that lam input is in nm, then converted to cm below.
      implicit none
      integer npass
      real T,l,h,c,k,temp(2),Pi,lam(npass),pass(npass),i,dl
      data k /1.380658e-16/
      data c /2.99792458e10/
      data h /6.6260755e-27/
      data Pi /3.141592654/

      planck=0. !initialize for integration
      do 10 i=2,npass
       
        l=lam(i)*1.0e-7 !convert nm -> cm
        dl=(lam(i)-lam(i-1))*1.0e-7 !calculate delta lambda
      
        temp(1)=2*Pi*h*c*c/(l**5) !planck at specific wavelength 
        temp(2)=exp(h*c/(l*k*T))-1.0
      
        planck=planck+temp(1)/temp(2)*dl*(pass(i)+pass(i-1))/2.0
 10   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readKeplerbandpass(nKepler,Keplerlam,Keplerpass)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nKepler,i,tempi
      real Keplerlam(nKepler),Keplerpass(nKepler)
      character*80 filename
      
      filename='/home/rowe/p555/ptempmod/kepler_trans.dat'
      
      open(unit=11,file=filename,status='old',err=901)
      
      do 10 i=1,nKepler
        read(11,*) tempi,Keplerpass(i)
        Keplerlam(i)=real(tempi)/10.0 !convert from A to nm
 10   continue
 
      close(11)
      goto 999
 901  write(0,*) "Cannot open: ",filename
      stop
      goto 999
 999  return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER n
      real adev,ave,curt,sdev,skew,var,data(n)
      INTEGER j
      real p,s,ep
      if(n.le.1)pause 'n must be at least 2 in moment'
      s=0.
      do 11 j=1,n
         s=s+data(j)
 11   enddo
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         adev=adev+abs(s)
         p=s*s
         var=var+p
         p=p*s
         skew=skew+p
         p=p*s
         curt=curt+p
 12   enddo
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=sqrt(var)
      if(var.ne.0.)then
         skew=skew/(n*sdev**3)
         curt=curt/(n*var**2)-3.
      else
         pause 'no skew or kurtosis when zero variance in moment'
      endif
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function ranb(Per,rhostar,seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer seed
      real Per,rhostar,b,Psec,adrs,G,Pi,cincl,ran2,incl,rnd
      
      Pi=acos(-1.e0)!define Pi and 2*Pi
      G=6.674e-11 !N m^2 kg^-2  Gravitation constant
      b=1.0
      
      Psec=Per*24.0*60.0*60.0 !sec ; period of planet
      adrs=1000.0*rhostar*G*(Per*86400.0e0)**2/(3.0e0*Pi)
      adrs=adrs**(1.0e0/3.0e0)
      cincl=b/adrs !cos(i)
      
      rnd=ran2(seed)
      incl=(pi/2.0-acos(cincl))*rnd
      ranb=cos(pi/2.0-incl)*adrs
      
c      write(0,*) ranb,incl,cincl
c      read(5,*)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION ran2(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     * IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     * IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function transitdur(Per,rhostar,rdr,b) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real temp(4),Psec,Pi,adrs,rdr,cincl,Per,rhostar,b,G,bb
      
      Pi=acos(-1.d0)!define Pi and 2*Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      
      Psec=Per*24.0*60.0*60.0 !sec ; period of planet
      adrs=1000.0*rhostar*G*(Per*86400.0d0)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0)
      cincl=b/adrs !cos(i)
      bb=b*b
      
c       write(0,*) Psec,adrs,cincl
        temp(1)=Psec/Pi
        temp(2)=1.0d0/adrs
        temp(3)=(1+rdr)**2.0-bb
        temp(4)=1-cincl*cincl
C     Transit duration in hours
        transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))/3600.0
 
      return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine histogram(nbin,npt,rp,nbinmax,bdatax,bdatay,title)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,i,nbinmax
      real rp(npt),datamin,datamax,bmax,bdatax(nbinmax),bdatay(nbinmax)
      character*80 title
      
C     Find datarange
      datamin=rp(1)
      datamax=rp(1)
      do 12 i=2,npt
        datamin=min(rp(i),datamin)
        datamax=max(rp(i),datamax)
 12   continue
 
      call bindata(nbin,npt,rp,bdatax,bdatay,datamin,datamax,bmax)
      
      call pgpage() !fresh plotting surface
      call pgslw(3)
      call pgsch(2.0)
c         call windowsetup(xb1,xb2,yb1,yb2) !make a square plotting surface
c         call pgvport(xb1,xb2,yb1,yb2)
c      call pgvport(0.2,1.0,0.2,0.9)
c      fc=2.*(datamax-datamin)/real(nbin) !center histograms
      call pgwindow(datamin,datamax,0.,bmax+0.1*bmax) !set size
         
C        Add axis labels
      call pglabel(title,"Relative Probability","")
      call pgbin(nbin,bdatax,bdatay,.false.) !plot the histogram
      
C     Shift axis scale to account for average removal
      call pgwindow(datamin,datamax,0.,1.0+0.1*1.0)
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0) !add boarders
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindata(nbin,npt,pdata,bdatax,bdatay,datamin,datamax,
     .   bmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nbin,npt,i,bin
      real pdata(npt),bdatax(nbin),bdatay(nbin),datamin,datamax,
     .   binsize,bmax,tsum

C     get size of bin.
      binsize=(datamax-datamin)/real(nbin-1)
      
C     initalize bdata
      do 20 i=1,nbin
         bdatay(i)=0.
 20   continue

C     loop over data and place data into proper bin.
      do 10 i=1,npt
C        calculate bin number
         bin=int((pdata(i)-datamin)/binsize)+1
         if((bin.gt.0).and.(bin.le.nbin)) bdatay(bin)=bdatay(bin)+1.0
 10   continue

C     get the max value in a bin and assign bin value.
      tsum=0.
      bmax=0.
      do 30 i=1,nbin
         bmax=max(bdatay(i),bmax)
         bdatax(i)=datamin+real(i-1)*binsize !added -1
         tsum=tsum+bdatay(i)
 30   continue
c      write(6,*) "Sum:",tsum
 
      return
      end