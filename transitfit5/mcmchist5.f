      program mcmchist5
      implicit none
      integer iargc,nmax,np,nunit,nbin,nbinmax,i,nfitm,nplanet,nfit,ii,
     .  npt,j,k,seed,now(3),jj
      parameter(nmax=5000000,nbinmax=500,nfitm=112)
      integer nd(nmax)
      real rp(nmax),bdatax(nbinmax),bdatay(nbinmax),errs(6),ave,std
      double precision mstar(nmax),age(nmax),z(nmax),rstar(nmax),
     .  rhostar(nmax),temp(nmax),lum(nmax),work(nmax),sol(nfitm),
     .  serr(nfitm,2),err(nfitm),dd(nmax),dumr,ran2,med,output(17)
      character*80 rhofile,title,parsfile,titles(18),mcmcfile,cline,
     .  names(18),name

      call itime(now)
      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      dumr=ran2(-seed)

      titles(1)="Stellar Density (g/cm\u3\d)"
      titles(2)='Limbdarkening(1)'
      titles(3)='Limbdarkening(2)'
      titles(4)='Limbdarkening(3)'
      titles(5)='Limbdarkening(4)'
      titles(6)='Dilution'
      titles(7)="Velocity Offset (m/s)"
      titles(8)='Phot zpt (ppm)'
      titles(9)='Epoch (BJD-2454900)'
      titles(10)='P-P0 (10\u-6\d days)'
      titles(11)='b'
      titles(12)="Rp/R*"
      titles(13)="e cos w"!"Eccentricity"
      titles(14)="e sin w"!"Argument of Pericenter (deg)"
      titles(15)="K (m/s)"
      titles(16)='Secondary Eclipse Depth (ppm)'
      titles(17)='Ellipsodial Amplitude (ppm)'
      titles(18)='Phase changes (ppm)'

      names(1)="rhostar (g/cm^3)"
      names(2)='NL1'
      names(3)='NL2'
      names(4)='NL3'
      names(5)='NL4'
      names(6)='Dilution'
      names(7)='gamma (m/s)'
      names(8)='Phot zpt (ppm)'
      names(9)='T0 (BJD-2454900)'
      names(10)='Period (days)'
      names(11)='b'
      names(12)="Rp/R*"
      names(13)="e cos w"
      names(14)="e sin w"
      names(15)="K (m/s)"
      names(16)='Occultation (ppm)'
      names(17)='Ellipsodial (ppm)'
      names(18)='Phases (ppm)'

      if(iargc().lt.2) goto 901 !check number of input parameters

      nbin=30 !change to optional input parameter and check nbinmax
      if(iargc().ge.3)then
        call getarg(3,cline)
        read(cline,*) nbin
      endif
      if(nbin.lt.2) goto 905
      if(nbin.gt.nbinmax) goto 906

      call pgopen('?') !open PGPlot device
c      call pgopen('/null')
      call pgask(.true.) !don't ask for new page.. just do it.
      call PGPAP ( 8.0 ,1.0) !paper size
      call pgsubp(3,3)  !break up plot into grid

      write(6,501) "Parameter    ","Median    ", "Stdev     ",
     .  "+1 sig     ", "-1 sig     ", "+2 sig     ", "-2 sig     ",
     .  "+3 sig     ", "-3 sig     "
      write(6,502) "----------------------------------------------------
     .------------------------------------------------------------------
     .---"
 501  format(A18,8(A13))
 502  format(A122)

C     Read in best-fit parameter file as well.
      call getarg(1,parsfile) !get mcmc parameters from rhostar analysis
      nunit=10
      open(unit=nunit,file=parsfile,status='old',err=903)
      call getfitpars(nunit,nfitm,nplanet,sol,serr,err)
      close(nunit)
      nfit=nplanet*10+8
      write(0,*) "nPlanet: ",nplanet


C     Now we read in MCMC results, cycling through as needed.
      call getarg(2,mcmcfile) !get mcmc parameters from rhostar analysis
      nunit=10
      open(unit=nunit,file=mcmcfile,status='old',err=904)

      do 10 i=1,8
c        write(6,*) titles(i),serr(i,2)
         if(serr(i,2).ne.0.0d0)then
            call getdata(nunit,i,npt,dd)
            title=titles(i)
            name=names(i)
            if(i.eq.8)then
                do 14 j=1,npt
                    dd(j)=dd(j)*1.0d6
 14             continue
            endif
            call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,bdatay,
     .          title,ave,std,errs)
c            write(6,500) sol(i),ave,std,(errs(j),j=1,6)
            call writetable(ave,std,errs,name)
            if(i.eq.1)then
               output(1)=ave
               output(2)=errs(1)
               output(3)=errs(2)
            endif
            rewind(nunit)
         endif
 10   continue

      do 11 i=1,nplanet
        do 12 j=1,10
            ii=8+(i-1)*10+j
            if(serr(ii,2).ne.0.0d0)then
                title=titles(8+j)
                name=names(8+j)
                call getdata(nunit,ii,npt,dd)
                if(j.eq.2)then !for Period we need to remove the median
                    call rqsort(npt,dd,nd) !changed npt to k
                    jj=npt/2
                    if(jj.le.0) jj=1
                    med=dd(nd(jj))
                    do 15 jj=1,npt
                        dd(jj)=(dd(jj)-med)!*1.0d6
 15                 continue
                endif

                call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .              bdatay,title,ave,std,errs)
c                write(6,500) sol(ii),ave,std,(errs(k),k=1,6)
                if(j.le.2) then
                    call writetable(real(sol(ii)),std,errs,name)
                else
                    call writetable(ave,std,errs,name)
                endif
                if(j.eq.1)then
                  output(4)=ave
                  output(5)=std
                elseif(j.eq.2)then
                  output(6)=sol(ii)
                  output(7)=std
                elseif(j.eq.3)then
                  output(8)=ave
                  output(9)=errs(1)
                  output(10)=errs(2)
                elseif(j.eq.4)then
                  output(11)=ave
                  output(12)=errs(1)
                  output(13)=errs(2)
                endif

                rewind(nunit)
            endif
 12     continue


        call getadrs(nunit,i,npt,dd,np,rstar,mstar,rhostar,seed)
        title="a/R\d\(2281)\u"
        call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .      bdatay,title,ave,std,errs)
c        write(6,500) ave,std,(errs(k),k=1,6)
        name="a/R*"
        call writetable(ave,std,errs,name)
        rewind(nunit)

c        call getasemi(nunit,i,npt,dd,np,mstar,seed)
c        title="a (AU)"
c        call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
c     .      bdatay,title,ave,std,errs)
cc        write(6,500) ave,std,(errs(k),k=1,6)
c        name="a (AU)"
c        call writetable(ave,std,errs,name)
c        rewind(nunit)

        call gettdep(nunit,i,npt,dd)
        title="Transit Depth (ppm)"
        call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .      bdatay,title,ave,std,errs)
c        write(6,500) ave,std,(errs(k),k=1,6)
        name="Tdepth (ppm)"
        call writetable(ave,std,errs,name)
        rewind(nunit)
        output(14)=ave
        output(15)=std

c        if((serr(8+10*(i-1)+7,2).ne.0.0d0).and.
c     .     (serr(8+10*(i-1)+4,2).ne.0.0d0))then
c            call getrhoplanet(nunit,i,npt,dd,np,mstar,rstar,seed)
c            title="Planet Density (g/cm\u3\d)"
c            call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
c     .          bdatay,title,ave,std,errs)
cc            write(6,500) ave,std,(errs(k),k=1,6)
c            name="rho_p (g/cm^3)"
c            call writetable(ave,std,errs,name)
c            rewind(nunit)
c        endif

c        call gett12(nunit,i,nmax,npt,dd,np,rstar,mstar,seed)
c        title="T\d12\u (h)"
c        call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
c     .      bdatay,title,ave,std,errs)
cc        write(6,500) ave,std,(errs(k),k=1,6)
c        name="T12 (h)"
c        call writetable(ave,std,errs,name)
c        rewind(nunit)

        call gettdur2(nunit,i,nmax,npt,dd,np,rstar,mstar,seed)
        title="T\dd\u (h)"
        call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .      bdatay,title,ave,std,errs)
c        write(6,500) ave,std,(errs(k),k=1,6)
        name="Tdur (h)"
        call writetable(ave,std,errs,name)
        rewind(nunit)
        output(16)=ave
        output(17)=std

 11   continue

 25   continue

 900  close(nunit)
      call pgclos()

      write(6,503) (output(i),i=1,17)
 503  format(17(1PE17.10,1X))

      goto 999
 901  write(0,*) "Usage: mcmchist5 n1.dat mcmc.dat"
      goto 999
 902  write(0,*) "Cannot open ",rhofile
      goto 999
 903  write(0,*) "Cannot open ",parsfile
      goto 999
 904  write(0,*) "Cannot open ",mcmcfile
      goto 999
 905  write(0,*) "Error: Nbin must be at least 2"
      goto 999
 906  write(0,*) "Increase nbinmax to at least",nbin
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine writetable(ave,std,errs,name)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i
      real ave,std,errs(6),minerr
      character*80 name

c      digit=log10(abs(errs(1)))
c      digit=min(log10(abs(errs(2))),digit)
c      ndigit=int(digit)
c      ndigit=ndigit-2
c      if(ndigit.lt.1) ndigit=ndigit-2
c      if(ndigit.ge.1) ndigit=ndigit-1
c      write(0,*) errs(1),errs(2)
c      write(0,*) ndigit,digit,log10(abs(errs(1))),log10(abs(errs(2)))
c      write(0,*) "ndigit:",ndigit

      minerr=min(errs(1),-errs(2))

      if(minerr.le.1.0e-5)then
        write(6,500) name,ave,std,(errs(i),i=1,6)
 500    format(A18,8(F12.7,1X))
      elseif(minerr.le.1.0e-4)then
        write(6,501) name,ave,std,(errs(i),i=1,6)
 501    format(A18,8(F11.6,2X))
      elseif(minerr.le.1.0e-3)then
        write(6,502) name,ave,std,(errs(i),i=1,6)
 502    format(A18,8(F10.5,3X))
      elseif(minerr.le.1.0e-2)then
        write(6,503) name,ave,std,(errs(i),i=1,6)
 503    format(A18,8(F9.4,4X))
      elseif(minerr.le.1.0e-1)then
        write(6,504) name,ave,std,(errs(i),i=1,6)
 504    format(A18,8(F8.3,5X))
      elseif(minerr.le.1.0)then
        write(6,505) name,ave,std,(errs(i),i=1,6)
 505    format(A18,8(F7.2,6X))
      elseif(minerr.le.10.0)then
        write(6,506) name,ave,std,(errs(i),i=1,6)
 506    format(A18,8(F6.1,7X))
      else
        if(ave.lt.10000.0)then
            write(6,507) name,ave,std,(errs(i),i=1,6)
 507        format(A18,8(F5.0,8X))
        else
            write(6,508) name,ave,std,(errs(i),i=1,6)
 508        format(A18,8(F7.0,6X))
        endif
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gettdep(nunit,nplanet,npt,tdep)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,j,col,i,dtype(1),nfitm,nplanetmax,nfit,
     .   nmax
      parameter(nfitm=112,nplanetmax=10,nmax=1500000)
      integer ntt(nplanetmax)
      double precision tobs(nplanetmax,nmax),omc(nplanetmax,nmax)
      double precision tdepth(npt),sol(nfitm),dumr,tmodel(1),tdep(npt),
     .   itime(1),time(1)

      nfit=nfitm
      itime(1)=1765.5/86400.0d0
      dtype(1)=0
      ntt(1)=0

C     Get all the parameters for a transit model
      col=10*(nplanet-1)
      i=1 !counter
      read(10,*) dumr
 10   read(nunit,*,end=11) dumr,dumr,dumr,(sol(j),j=1,8),
     .   (dumr,j=1,col),(sol(j),j=9,18)

        time(1)=sol(9)

C     FIX ME.... HERE
c        call transitmodel(nfit,1,sol,1,time,itime,tmodel,dtype)
        call transitmodel(nfit,1,nplanetmax,sol,nmax,1,time,
     .  itime,ntt,tobs,omc,tmodel,dtype)

        tdep(i)=(1.0d0-tmodel(1)+sol(8))*1.0d6
c        write(0,*) (1.0d0-tmodel)*1.0d6
c        read(5,*)

        i=i+1
        goto 10
 11   continue
      npt=i-1 !update counter

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getrhoplanet(nunit,nplanet,npt,mp,np,mstar,rstar,seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,np,seed,i,j,k,col
      double precision mp(npt),rstar(np),mstar(np),ecw,esw,per,bb,Kr,b,
     .  Pi,Msun,Rsun,G,aConst,Psec,M1,R1,incl,dumr,ran2,asemi,ac,bc,C1,
     .  C2,ct(3),Mearth,rdr,pmass,prad,fourthirdsPi

      ct(1)=1.0d0/3.0d0
      ct(2)=2.0d0**ct(1)
      ct(3)=3.0d0*sqrt(3.0d0)

      Pi=acos(-1.d0)   !Pi
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.9742d24 !kg mass of Earth
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      fourthirdsPi=4.0d0/3.0d0*Pi

C     Need period from transit models
      col=8+10*(nplanet-1)+2
      i=1 !counter
      read(10,*) dumr
 10   read(nunit,*,end=11) (dumr,j=1,col+2),per,b,rdr,ecw,esw,Kr

c        Kr=Kr/1.9103

        bb=b*b !b to b^2
        Psec=per*8.64d4 !sec ; period of planet
        k=ran2(seed)*(np-1)+1 !pick a random selection from isochrones
        M1=mstar(k)*Msun
        R1=rstar(k)*Rsun
        asemi=(M1)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst
        incl=acos(b*R1/asemi) !inclination in radians

        bc=Kr**3*Psec/(2*Pi*G*sin(incl)**3)
        ac=M1

        if(bc.gt.0)then
            C2=sqrt(27.0*ac**4*bc**2+4.0d0*ac**3*bc**3)
            C1=(27.0*ac**2*bc+ct(3)*C2+18.0d0*ac*bc**2+2.0d0*bc**3)**
     .          ct(1)
            pmass=C1/(3.0d0*ct(2))-ct(2)*(-6.0d0*ac*bc-bc**2)/
     .          (3.0d0*C1)+bc/3.0d0
c            mp(i)=mp(i)/Mearth !convert to Earth Units
        else
            goto 10 !we don't count negative mass
        endif

        prad=R1*rdr !planet radius

        mp(i)=pmass/(fourthirdsPi*prad*prad*prad)/1000.0d0

c        if((mp(i).gt.-9.9e30).and.(mp(i).lt.9.9e30))then
c            continue
c        else
c            write(0,500) b,Psec,Kr,M1,R1,asemi,sin(incl)
c            write(0,500) ac,bc,C2,C1
c            write(0,500) C1/(3.0d0*ct(2)),
c     .                   ct(2)*(-6.0d0*ac*bc-bc**2)/3.0d0*C1,
c     .                   bc/3.0d0
c            write(0,500) mp(i)
c            read(5,*)
c 500        format(10(1PE16.9,1X))
c        endif

        i=i+1
        goto 10
 11   continue
      npt=i-1 !update counter

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getpmass(nunit,nplanet,npt,mp,np,rstar,mstar,seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,np,seed,i,j,k,col
      double precision mp(npt),rstar(np),mstar(np),ecw,esw,per,bb,Kr,b,
     .  Pi,Msun,Rsun,G,aConst,Psec,M1,R1,incl,dumr,ran2,asemi,ac,bc,C1,
     .  C2,ct(3),Mearth

      ct(1)=1.0d0/3.0d0
      ct(2)=2.0d0**ct(1)
      ct(3)=3.0d0*sqrt(3.0d0)

      Pi=acos(-1.d0)   !Pi
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.9742d24 !kg mass of Earth
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

C     Need period from transit models
      col=8+10*(nplanet-1)+2
      i=1 !counter
      read(10,*) dumr
 10   read(nunit,*,end=11) (dumr,j=1,col+2),per,b,dumr,ecw,esw,Kr

c        Kr=Kr/1.9103

        bb=b*b
        Psec=per*8.64d4 !sec ; period of planet
        k=ran2(seed)*(np-1)+1 !pick a random selection from isochrones
        M1=mstar(k)*Msun
        R1=rstar(k)*Rsun
        asemi=(M1)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst
        incl=acos(b*R1/asemi) !inclination in radians

        bc=Kr**3*Psec/(2*Pi*G*sin(incl)**3)
        ac=M1

        if(bc.gt.0)then
            C2=sqrt(27.0*ac**4*bc**2+4.0d0*ac**3*bc**3)
            C1=(27.0*ac**2*bc+ct(3)*C2+18.0d0*ac*bc**2+2.0d0*bc**3)**
     .          ct(1)
            mp(i)=C1/(3.0d0*ct(2))-ct(2)*(-6.0d0*ac*bc-bc**2)/
     .          (3.0d0*C1)+bc/3.0d0
            mp(i)=mp(i)/Mearth !convert to Earth Units
        else
            goto 10 !we don't count negative mass
        endif

c        if((mp(i).gt.-9.9e30).and.(mp(i).lt.9.9e30))then
c            continue
c        else
c            write(0,500) b,Psec,Kr,M1,R1,asemi,sin(incl)
c            write(0,500) ac,bc,C2,C1
c            write(0,500) C1/(3.0d0*ct(2)),
c     .                   ct(2)*(-6.0d0*ac*bc-bc**2)/3.0d0*C1,
c     .                   bc/3.0d0
c            write(0,500) mp(i)
c            read(5,*)
c 500        format(10(1PE16.9,1X))
c        endif

        i=i+1
        goto 10
 11   continue
      npt=i-1 !update counter

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gett12(nunit,nplanet,nmax,npt,t12,np,rstar,mstar,
     .  seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,np,seed,i,j,k,col,nmax
      double precision tdur,rstar(np),mstar(np),Psec,per,dumr,Pi,
     .  Msun,Rsun,G,aConst,temp(5),ran2,M1,R1,asemi,bb,b,cincl,rdr,R2,
     .  adrs,rhostar,tF,t12(nmax)

      Pi=acos(-1.d0)   !Pi
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

C     Need period from transit models
      col=8+10*(nplanet-1)+2
      i=1 !counter
      read(10,*) dumr
c 10   read(nunit,*,end=11) (dumr,j=1,col+2),per,bb,rdr
 10   read(nunit,*,end=11) dumr,dumr,dumr,rhostar,(dumr,j=1,col-2),per,
     .      b,rdr
        bb=b*b
        k=ran2(seed)*(np-1)+1 !pick a random selection from isochrones
        Psec=per*8.64d4 !sec ; period of planet
        adrs=1000.0*rhostar*G*(per*86400.0d0)**2/(3.0d0*Pi)
        adrs=adrs**(1.0d0/3.0d0) !a/R*
        cincl=b/adrs !cos(i)

        temp(1)=Psec/Pi
        temp(2)=1.0d0/adrs
        temp(3)=(1+rdr)**2.0-bb
        temp(4)=1-cincl*cincl
        temp(5)=(1-rdr)**2.0-bb
C     Transit duration in hours
        tdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))/3600.0
        tF=temp(1)*asin(temp(2)*sqrt(temp(5)/temp(4)))/3600.0
        t12(i)=(tdur-tF)/2.0d0
        if((t12(i).gt.-9.9e30).and.(t12(i).lt.9.9e30))then
c            write(0,500) per,bb,rdr,tdur(i)
            i=i+1
        else
            write(0,500) per,bb,rdr
            write(0,500) Psec,adrs,cincl
            write(0,500) (temp(j),j=1,4)
            write(0,500) t12(i),dble(i)
c            read(5,*)
 500        format(7(1PE16.9,1X))
        endif

        goto 10
 11   continue
      npt=i-1

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gettdur2(nunit,nplanet,nmax,npt,tdur,np,rstar,mstar,
     .  seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,np,seed,i,j,k,col,nmax
      double precision tdur(nmax),rstar(np),mstar(np),Psec,per,dumr,Pi,
     .  Msun,Rsun,G,aConst,temp(4),ran2,M1,R1,asemi,bb,b,cincl,rdr,R2,
     .  adrs,rhostar

      Pi=acos(-1.d0)   !Pi
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

C     Need period from transit models
      col=8+10*(nplanet-1)+2
      i=1 !counter
      read(10,*) dumr
c 10   read(nunit,*,end=11) (dumr,j=1,col+2),per,bb,rdr
 10   read(nunit,*,end=11) dumr,dumr,dumr,rhostar,(dumr,j=1,col-2),per,
     .      b,rdr
        bb=b*b
        k=ran2(seed)*(np-1)+1 !pick a random selection from isochrones
        Psec=per*8.64d4 !sec ; period of planet
        adrs=1000.0*rhostar*G*(per*86400.0d0)**2/(3.0d0*Pi)
        adrs=adrs**(1.0d0/3.0d0) !a/R*
        cincl=b/adrs !cos(i)

        temp(1)=Psec/Pi
        temp(2)=1.0d0/adrs
        temp(3)=(1+rdr)**2.0-bb
        temp(4)=1-cincl*cincl
C     Transit duration in hours
        tdur(i)=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))/3600.0
        if((tdur(i).gt.-9.9e30).and.(tdur(i).lt.9.9e30))then
c            write(0,500) per,bb,rdr,tdur(i)
            i=i+1
        else
            write(0,500) per,bb,rdr
            write(0,500) Psec,adrs,cincl
            write(0,500) (temp(j),j=1,4)
            write(0,500) tdur(i),dble(i)
c            read(5,*)
 500        format(7(1PE16.9,1X))
        endif

        goto 10
 11   continue
      npt=i-1

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gettdur(nunit,nplanet,nmax,npt,tdur,np,rstar,mstar,
     .  seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,np,seed,i,j,k,col,nmax
      double precision tdur(nmax),rstar(np),mstar(np),Psec,per,dumr,Pi,
     .  Msun,Rsun,G,aConst,temp(4),ran2,M1,R1,asemi,bb,b,cincl,rdr,R2

      Pi=acos(-1.d0)   !Pi
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

C     Need period from transit models
      col=8+10*(nplanet-1)+2
      i=1 !counter
      read(10,*) dumr
 10   read(nunit,*,end=11) (dumr,j=1,col+2),per,b,rdr
        bb=b*b
        k=ran2(seed)*(np-1)+1 !pick a random selection from isochrones
        Psec=per*8.64d4 !sec ; period of planet
        M1=mstar(k)*Msun
        R1=rstar(k)*Rsun
        R2=R1*rdr
        asemi=(M1)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst
        cincl=b*R1/asemi !cos(i)

        temp(1)=Psec/Pi
        temp(2)=R1/asemi
        temp(3)=(1+(R2/R1))**2.0-((asemi/R1)*cincl)**2.0
        temp(4)=1-cincl*cincl
C     Transit duration in hours
        tdur(i)=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))/3600.0
        if((tdur(i).gt.-9.9e30).and.(tdur(i).lt.9.9e30))then
c            write(0,500) per,bb,rdr,tdur(i)
            i=i+1
        else
            write(0,500) per,bb,rdr
            write(0,500) Psec,M1,R1,R2,asemi,cincl,rdr
            write(0,500) (temp(j),j=1,4)
            write(0,500) tdur(i),dble(i)
c            read(5,*)
 500        format(7(1PE16.9,1X))
        endif

        goto 10
 11   continue
      npt=i-1

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getasemi(nunit,nplanet,npt,asemi,np,mstar,seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,npt,np,seed,nplanet,i,j,k,col
      double precision asemi(npt),mstar(np),per,Pi,Msun,
     .  G,aConst,M1,Psec,dumr,ran2,AU

      Pi=acos(-1.d0)   !Pi
      Msun=1.9891d30 !kg  mass of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      AU=1.49598e11

C     Need period from transit models
      col=8+10*(nplanet-1)+2
      i=1 !counter
      read(10,*) dumr
 10   read(nunit,*,end=11) (dumr,j=1,col+2),per
        k=ran2(seed)*(np-1)+1 !pick a random selection from isochrones
        Psec=per*8.64d4 !sec ; period of planet
        M1=mstar(k)*Msun
        asemi(i)=(M1)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst/AU
        i=i+1
        goto 10
 11   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getadrs(nunit,nplanet,npt,adrs,np,rstar,mstar,rhostar,
     .  seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,npt,np,seed,nplanet,i,j,k,col
      double precision adrs(npt),rstar(np),mstar(np),per,Pi,Msun,Rsun,
     .  G,aConst,M1,R1,Psec,asemi,dumr,ran2,rhostar(np),rhostarmodel

      Pi=acos(-1.d0)   !Pi
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

C     Need period from transit models
      col=8+10*(nplanet-1)+2
      i=1 !counter
      read(10,*) dumr
 10   read(nunit,*,end=11) dumr,dumr,dumr,rhostarmodel,(dumr,j=1,col-2),
     .  per
         k=ran2(seed)*(np-1)+1 !pick a random selection from isochrones
c        Psec=per*8.64d4 !sec ; period of planet
c        M1=mstar(k)*Msun
c        R1=rstar(k)*Rsun
c        asemi=(M1)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst
c        adrs(i)=asemi/R1


c        adrs(i)=1000.0*rhostar(k)*G*(per*86400.0d0)**2/(3.0d0*Pi)
        adrs(i)=1000.0*rhostarmodel*G*(per*86400.0d0)**2/(3.0d0*Pi)
        adrs(i)=adrs(i)**(1.0d0/3.0d0)
c        write(0,*) rhostar,per
c        read(5,*)

        i=i+1
        goto 10
 11   continue
      npt=i-1 !update counter

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getincl(nunit,nplanet,npt,incl,np,rstar,mstar,seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,seed,np,col,i,j,k
      double precision incl(npt),rstar(np),mstar(np),ran2,per,bb,b,Pi,
     .  Msun,Rsun,G,aConst,asemi,M1,Psec,R1,dumr

      Pi=acos(-1.d0)   !Pi
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

C     Need impact parameter, period from transit models
      col=8+10*(nplanet-1)+2
      i=1 !counter
      read(10,*) dumr
 10   read(nunit,*,end=11) (dumr,j=1,col+2),per,b
        k=ran2(seed)*(np-1)+1 !pick a random selection from isochrones
        bb=b*b
        M1=mstar(k)*Msun
        R1=rstar(k)*Rsun
        Psec=per*8.64d4 !sec ; period of planet
        asemi=(M1)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst
        incl(i)=acos(b*R1/asemi)*180.0/Pi !inclination in degrees
        i=i+1
        goto 10
 11   continue
      npt=i-1 !update counter

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getprad(nunit,nplanet,npt,rp,np,rstar,seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,np,nplanet,col,i,j,k,seed,nunit
      double precision rp(npt),rstar(np),dumr,rdr,ran2

C     Read in r/R* column
      col=8+10*(nplanet-1)+4
      i=1 !counter
      read(10,*) dumr
 10   read(nunit,*,end=11) (dumr,j=1,col+2),rdr
        k=ran2(seed)*(np-1)+1 !pick a random selection from rstar
        rp(i)=rstar(k)*rdr*109.17
        i=i+1
        goto 10
 11   continue
      npt=i-1

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getlogg(np,logg,mstar,rstar)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer np,i
      double precision logg(np),mstar(np),rstar(np),Msun,Rsun,G

      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant

      do 10 i=1,np
        logg(i)=log10(1.0e2*G*Msun*mstar(i)/
     .      (Rsun*Rsun*rstar(i)*rstar(i)))
  10  continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getdata(nunit,col,npt,dd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer col,npt,i,j,k,nunit
      double precision dd(npt),dumr

      do 5 i=1,1!5000
      read(nunit,*) dumr
 5    continue

      k=0
      i=1
  9   continue
 10   read(nunit,*,end=11) (dumr,j=1,col+2),dd(i)
c         write(6,*) dd(i)
c         if((col.eq.6).and.(dd(i).gt.89.95)) goto 10
c         if(col.eq.8) dd(i)=dd(i)*1.0e6

c         k=k+1
c         if(k.lt.10) goto 10
c         k=0

         i=i+1
         goto 10
 11   continue

      npt=i-1

 500  format(200(E14.7,1X))
      goto 999
 901  write(6,*) "Error on line: ",i
      goto 999

 999  return
      end

