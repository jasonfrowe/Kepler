      program mcmchist5
      implicit none
      integer iargc,nmax,np,nunit,nbin,nbinmax,i,nfitm,nplanet,nfit,ii,
     .  npt,j,k,seed,now(3),jj
      parameter(nmax=5000000,nbinmax=500,nfitm=112)
      integer nd(nmax)
      real rp(nmax),bdatax(nbinmax),bdatay(nbinmax),errs(6),ave,std
      double precision mstar(nmax),age(nmax),z(nmax),rstar(nmax),
     .  rhostar(nmax),temp(nmax),lum(nmax),work(nmax),sol(nfitm),
     .  serr(nfitm,2),err(nfitm),dd(nmax),dumr,ran2,med,output(17),
     .   outerr(17)
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
      
      if(iargc().lt.3) goto 901 !check number of input parameters
      
      nbin=30 !change to optional input parameter and check nbinmax
      if(iargc().ge.4)then
        call getarg(4,cline)
        read(cline,*) nbin
      endif
      if(nbin.lt.2) goto 905
      if(nbin.gt.nbinmax) goto 906
      
      call getarg(2,rhofile) !get mcmc parameters from rhostar analysis.
      nunit=10
      open(unit=nunit,file=rhofile,status='old',err=902)
      call getrhostar(nunit,np,nmax,mstar,age,z,rstar,rhostar,temp,lum)
      close(nunit)

      call pgopen('?') !open PGPlot device
      call pgask(.false.) !don't ask for new page.. just do it.
      call PGPAP ( 8.0 ,1.0) !paper size
      call pgsubp(5,5)  !break up plot into grid
      
      write(6,501) "Parameter    ","Median    ", "Stdev     ", 
     .  "+1 sig     ", "-1 sig     ", "+2 sig     ", "-2 sig     ", 
     .  "+3 sig     ", "-3 sig     "
      write(6,502) "----------------------------------------------------
     .------------------------------------------------------------------
     .---"
 501  format(A18,8(A13))   
 502  format(A122)  
            
C     Now lets make histograms of the stellar-parameters.
      title='Stellar Mass (M\d\(2281)\u)'
      call histogram(np,rp,mstar,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="M* (Msun)"
      call writetable(ave,std,errs,name)
  500 format(9(F11.6,1X))
        
      title='Stellar Radius (R\d\(2281)\u)'
      call histogram(np,rp,rstar,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="R* (Rsun)"
      call writetable(ave,std,errs,name)
     
      title='Z'
      call histogram(np,rp,z,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="Z"
      call writetable(ave,std,errs,name)
     
      title='T\deff\u (K)'
      call histogram(np,rp,temp,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="Teff (K)"
      call writetable(ave,std,errs,name)
            
      title='Stellar Luminosity (L\d\(2281)\u)'
      call histogram(np,rp,lum,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="L (Lsun)"
      call writetable(ave,std,errs,name)
      
      title='Age (Gyr)'
      call histogram(np,rp,age,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="Age (Gyr)"
      call writetable(ave,std,errs,name)
  
C     logg here
      call getlogg(np,dd,mstar,rstar)
      title='log(g)'
      call histogram(np,rp,dd,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="log(g) (cgs)"
      call writetable(ave,std,errs,name)

C     rhostar here
      title='Stellar Density (g/cm\u3\d)'
      call histogram(np,rp,rhostar,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6) 
      name="rhostar (g/cm^3)"
      call writetable(ave,std,errs,name)
      
C     Read in best-fit parameter file as well.
      call getarg(1,parsfile) !get mcmc parameters from rhostar analysis
      nunit=10
      open(unit=nunit,file=parsfile,status='old',err=903)
      call getfitpars(nunit,nfitm,nplanet,sol,serr,err)
      close(nunit)
      nfit=nplanet*10+8
      write(0,*) "nPlanet: ",nplanet
      
      
C     Now we read in MCMC results, cycling through as needed.      
      call getarg(3,mcmcfile) !get mcmc parameters from rhostar analysis
      nunit=10
      open(unit=nunit,file=mcmcfile,status='old',err=904)

c      goto 50

c        write(0,*) "Tdur"
c        call gettdur(nunit,1,nmax,npt,dd,np,rstar,mstar,seed)
c        title="T\dd\u (h)"
c        call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
c     .      bdatay,title,ave,std,errs)
c        write(6,500) ave,std,(errs(k),k=1,6)
c        rewind(nunit)
c        write(0,*) "Tdur"

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
            rewind(nunit)
            if(i.eq.1) then
               output(1)=ave
               outerr(1)=std
            endif
         endif
 10   continue
      
      do 11 i=1,nplanet
         do 17 j=2,17
            output(j)=0.0d0 !arrays for machine readable output
            outerr(j)=0.0d0
 17      continue

c         if(i.gt.8)then !this is more linking titles
c            ii=i-10*((i-9)/10)
c         else
c            ii=i
c         endif
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

C              Changed from b^2 to b, so these lines are not needed.
c                if(j.eq.3)then
c                    do 13 k=1,npt
c                        dd(k)=sqrt(dd(k)) !convert b^2 to b
c 13                 continue
c                endif
                
c                if(j.eq.7)then  !fix Kr for Solar stars
c                    do 16 k=1,npt
c                        dd(k)=dd(k)/1.9103
c 16                 continue
c                endif
 
                call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .              bdatay,title,ave,std,errs)
c                write(6,500) sol(ii),ave,std,(errs(k),k=1,6)
                if(j.le.2) then
                    call writetable(real(sol(ii)),std,errs,name)
                else
                    call writetable(ave,std,errs,name)
                endif
                rewind(nunit)
                if(j.eq.1) then
                  output(2)=ave
                  outerr(2)=std
                elseif(j.eq.2)then
                  output(3)=real(sol(ii))
                  outerr(3)=std
                elseif(j.eq.3)then
                  output(5)=ave
                  outerr(5)=std
                elseif(j.eq.4)then
                  output(4)=ave
                  outerr(4)=std
                elseif(j.eq.5)then
                  output(7)=ave
                  outerr(7)=std
                elseif(j.eq.6)then
                  output(8)=ave
                  outerr(8)=std
                elseif(j.eq.7)then
                  output(6)=ave
                  outerr(6)=std
                endif

            endif
 12     continue

        if(serr(8+10*(i-1)+7,2).ne.0.0d0)then
            call getpmass(nunit,i,npt,dd,np,mstar,rstar,seed)
            title="M\dp\u (Mearth)"
            call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .          bdatay,title,ave,std,errs)
c            write(6,500) ave,std,(errs(k),k=1,6)
            name="Mp (Mearth)"
            call writetable(ave,std,errs,name)
            write(0,*) ave,std
            rewind(nunit)
            output(9)=ave
            outerr(9)=std
        endif
 
        if(serr(8+10*(i-1)+4,2).ne.0.0d0)then
            call getprad(nunit,i,npt,dd,np,rstar,seed)
            title="R\dp\u (Rearth)"
            call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .          bdatay,title,ave,std,errs)
c            write(6,500) ave,std,(errs(k),k=1,6)
            name="Rp (Rearth)"
            call writetable(ave,std,errs,name)
            rewind(nunit)
            output(10)=ave
            outerr(10)=std
        endif

        if(serr(8+10*(i-1)+3,2).ne.0.0d0)then
            call getincl(nunit,i,npt,dd,np,rstar,mstar,seed)
            title="i (deg)"
            call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .          bdatay,title,ave,std,errs)
c            write(6,500) ave,std,(errs(k),k=1,6)
            name="i (deg)"
            call writetable(ave,std,errs,name)
            rewind(nunit)
            output(11)=ave
            outerr(11)=std
        endif
        
        call getadrs(nunit,i,npt,dd,np,rstar,mstar,rhostar,seed)
        title="a/R\d\(2281)\u"
        call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .      bdatay,title,ave,std,errs)
c        write(6,500) ave,std,(errs(k),k=1,6)
        name="a/R*"
        call writetable(ave,std,errs,name)
        rewind(nunit)
        output(12)=ave
        outerr(12)=std
        
        call getasemi(nunit,i,npt,dd,np,mstar,seed)
        title="a (AU)"
        call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .      bdatay,title,ave,std,errs)
c        write(6,500) ave,std,(errs(k),k=1,6)
        name="a (AU)"
        call writetable(ave,std,errs,name)
        rewind(nunit)
        output(13)=ave
        outerr(13)=std
        
        if(serr(8+10*(i-1)+4,2).ne.0.0d0)then
         call gettdep(nunit,i,npt,dd)
         title="Transit Depth (ppm)"
         call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .      bdatay,title,ave,std,errs)
c        write(6,500) ave,std,(errs(k),k=1,6)
         name="Tdepth (ppm)"
         call writetable(ave,std,errs,name)
         rewind(nunit)
         output(14)=ave
         outerr(14)=std
        endif
        
c        if((serr(8+10*(i-1)+5,2).ne.0.0d0).and.
c     .     (serr(8+10*(i-1)+6,2).ne.0.0d0))then
c         call geteccn(nunit,i,npt,dd)
c         title="e"
c         call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
c     .      bdatay,title,ave,std,errs)
cc        write(6,500) ave,std,(errs(k),k=1,6)
c         name="e"
c         call writetable(ave,std,errs,name)
c         rewind(nunit)
c
c         call getomega(nunit,i,npt,dd)
c         title="w (deg)"
c         call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
c     .      bdatay,title,ave,std,errs)
cc        write(6,500) ave,std,(errs(k),k=1,6)
c         name="w (deg)"
c         call writetable(ave,std,errs,name)
c         rewind(nunit)
c        endif


        if((serr(8+10*(i-1)+7,2).ne.0.0d0).and.
     .     (serr(8+10*(i-1)+4,2).ne.0.0d0))then
            call getrhoplanet(nunit,i,npt,dd,np,mstar,rstar,seed)
            title="Planet Density (g/cm\u3\d)"
            call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .          bdatay,title,ave,std,errs)
c            write(6,500) ave,std,(errs(k),k=1,6)
            name="rho_p (g/cm^3)"
            call writetable(ave,std,errs,name)
            rewind(nunit)
            output(15)=ave
            outerr(15)=std

            call getloggplanet(nunit,i,npt,dd,np,mstar,rstar,seed)
            title="Planet log(g) (cgs)"
            call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .          bdatay,title,ave,std,errs)
c            write(6,500) ave,std,(errs(k),k=1,6)
            name="log(g)_p (cgs)"
            call writetable(ave,std,errs,name)
            rewind(nunit)
            output(16)=ave
            outerr(16)=std

        endif

c        if(serr(8+10*(i-1)+4,2).ne.0.0d0)then
            call getteq(nunit,i,npt,dd,np,mstar,rstar,temp,seed)
            title="T\deq\u (K)"
            call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .          bdatay,title,ave,std,errs)
c            write(6,500) ave,std,(errs(k),k=1,6)
            name="Teq (K)"
            call writetable(ave,std,errs,name)
            rewind(nunit)
            output(17)=ave
            outerr(17)=std
c        endif
 
         call getsrad(nunit,i,npt,dd,np,mstar,rstar,temp,seed)
         title="S (Searth)"
         call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
     .    bdatay,title,ave,std,errs)
c            write(6,500) ave,std,(errs(k),k=1,6)
         name="S (Searth)"
         call writetable(ave,std,errs,name)
         rewind(nunit)
         output(17)=ave
         outerr(17)=std

c        call gett12(nunit,i,nmax,npt,dd,np,rstar,mstar,seed)
c        title="T\d12\u (h)"
c        call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
c     .      bdatay,title,ave,std,errs)
cc        write(6,500) ave,std,(errs(k),k=1,6)
c        name="T12 (h)"
c        call writetable(ave,std,errs,name)
c        rewind(nunit)        
 
c        call gettdur2(nunit,i,nmax,npt,dd,np,rstar,mstar,seed)
c        title="T\dd\u (h)"
c        call histogram(npt,rp,dd,work,nd,nbin,nbinmax,bdatax,
c     .      bdatay,title,ave,std,errs)
cc        write(6,500) ave,std,(errs(k),k=1,6)
c        name="Tdur (h)"
c        call writetable(ave,std,errs,name)
c        rewind(nunit)

      write(6,503) (output(j),j=1,17)
      write(6,503) (outerr(j),j=1,17)
 503  format(F9.6,1X,F11.6,1X,F13.8,1X,F8.6,1X,F6.4,1X,F7.2,1X,
     .   2(F6.3,1X),F7.2,1X,F6.3,1X,F6.3,1X,F8.3,1X,F10.6,1X,F8.2,1X,
     .   F5.2,1X,F5.3,1X,F5.0)

 11   continue

 25   continue

 900  close(nunit)
      call pgclos()
      
      goto 999
 901  write(0,*) "Usage: mcmchist5 n1.dat rhoboot.file mcmc.dat"
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
      subroutine geteccn(nunit,nplanet,npt,eccn)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,j,col,i,dtype(1),nfitm,nplanetmax,nfit,
     .   nmax
      parameter(nfitm=112,nplanetmax=10,nmax=2000000)
      double precision eccn(npt),sol(nfitm),esw,ecw,dumr

      nfit=nfitm

C     Get all the parameters for a transit model
      col=10*(nplanet-1)
      i=1 !counter
      read(10,*) dumr
 10   read(nunit,*,end=11) dumr,dumr,dumr,(sol(j),j=1,8),
     .   (dumr,j=1,col),(sol(j),j=9,18)

        ecw=sol(14)
        esw=sol(13)
        eccn(i)=sqrt(ecw*ecw+esw*esw) !eccentricity
        if(eccn(i).ge.1.0) eccn=0.99

        i=i+1
        goto 10
 11   continue
      npt=i-1 !update counter

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getomega(nunit,nplanet,npt,w)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,j,col,i,dtype(1),nfitm,nplanetmax,nfit,
     .   nmax
      parameter(nfitm=112,nplanetmax=10,nmax=2000000)
      double precision eccn,sol(nfitm),esw,ecw,w(npt),Pi,tPi,dumr

      nfit=nfitm
      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi

C     Get all the parameters for a transit model
      col=10*(nplanet-1)
      i=1 !counter
      read(10,*) dumr
 10   read(nunit,*,end=11) dumr,dumr,dumr,(sol(j),j=1,8),
     .   (dumr,j=1,col),(sol(j),j=9,18)

        ecw=sol(14)
        esw=sol(13)
        eccn=sqrt(ecw*ecw+esw*esw) !eccentricity
        if(eccn.ge.1.0) eccn=0.99
        if(eccn.eq.0.0d0)then
            w(i)=0.0d0
        else
            if(ecw.eq.0.0d0)then
                w(i)=Pi/2.0d0
            else
                w(i)=atan(esw/ecw)
            endif
            if((ecw.gt.0.0d0).and.(esw.lt.0.0d0))then
                w(i)=tPi+w(i)
            elseif((ecw.lt.0.0d0).and.(esw.ge.0.0d0))then
                w(i)=Pi+w(i)
            elseif((ecw.le.0.0d0).and.(esw.lt.0.0d0))then
                w(i)=Pi+w(i)
            endif
        endif

        i=i+1
        goto 10
 11   continue
      npt=i-1 !update counter

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gettdep(nunit,nplanet,npt,tdep)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,j,col,i,dtype(1),nfitm,nplanetmax,nfit,
     .   nmax
      parameter(nfitm=112,nplanetmax=10,nmax=2000000)
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
      subroutine getsrad(nunit,nplanet,npt,mp,np,mstar,rstar,temp,seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,np,seed,col,i,j,k
      double precision mp(npt),mstar(np),rstar(np),temp(np),dumr,per,b,
     .   rdr,ecw,esw,Kr,Pi,Msun,Mearth,Rsun,G,aConst,fourthirdsPi,Psec,
     .   M1,R1,Teff,asemi,ran2,Tsun,AU

      Pi=acos(-1.d0)   !Pi
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.9742d24 !kg mass of Earth
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      fourthirdsPi=4.0d0/3.0d0*Pi
      AU=1.49597871d11
      Tsun=5781.6d0

C     Need period from transit models
      col=8+10*(nplanet-1)+2
      i=1 !counter
      read(nunit,*) dumr
 10   read(nunit,*,end=11) (dumr,j=1,col+2),per,b,rdr,ecw,esw,Kr

         Psec=per*8.64d4 !sec ; period of planet
         k=ran2(seed)*(np-1)+1 !pick a random selection from isochrones
         M1=mstar(k)*Msun
         R1=rstar(k)*Rsun
         Teff=temp(k)
         asemi=(M1)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst

         mp(i)=rstar(k)*rstar(k)*(Teff/Tsun)**4.0*(AU/asemi)**2.0

        i=i+1
        goto 10
 11   continue
      npt=i-1 !update counter

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getteq(nunit,nplanet,npt,mp,np,mstar,rstar,temp,seed)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nplanet,npt,np,seed,i,j,k,col
      double precision mp(npt),rstar(np),mstar(np),ecw,esw,per,bb,Kr,b,
     .  Pi,Msun,Rsun,G,aConst,Psec,M1,R1,incl,dumr,ran2,asemi,ac,bc,C1,
     .  C2,ct(3),Mearth,rdr,pmass,prad,fourthirdsPi,Teff,temp(np),f,Ab


c      write(0,*) "hello1"

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
      read(nunit,*) dumr
 10   read(nunit,*,end=11) (dumr,j=1,col+2),per,b,rdr,ecw,esw,Kr

c        Kr=Kr/1.9103

        bb=b*b !b to b^2
        Psec=per*8.64d4 !sec ; period of planet
        k=ran2(seed)*(np-1)+1 !pick a random selection from isochrones
        M1=mstar(k)*Msun
        R1=rstar(k)*Rsun
        Teff=temp(k)
        asemi=(M1)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst
        incl=acos(b*R1/asemi) !inclination in radians

        bc=Kr**3*Psec/(2*Pi*G*sin(incl)**3)
        ac=M1

c        if(bc.gt.0)then
c            C2=sqrt(27.0*ac**4*bc**2+4.0d0*ac**3*bc**3)
c            C1=(27.0*ac**2*bc+ct(3)*C2+18.0d0*ac*bc**2+2.0d0*bc**3)**
c     .          ct(1)
c            pmass=C1/(3.0d0*ct(2))-ct(2)*(-6.0d0*ac*bc-bc**2)/
c     .          (3.0d0*C1)+bc/3.0d0
cc            mp(i)=mp(i)/Mearth !convert to Earth Units
c        else
c            goto 10 !we don't count negative mass
c        endif

c        prad=R1*rdr !planet radius

c        mp(i)=pmass/(fourthirdsPi*prad*prad*prad)/1000.0d0
c         mp(i)=log10(G*pmass/(prad*prad))+2.0d0 !add 2 to make cgs
         f=ran2(seed)+1.0d0
         Ab=ran2(seed)/2.0d0
         mp(i)=Teff*sqrt(R1/(2.0d0*asemi))*(f*(1.0-Ab))**0.25
c         write(0,*) mp(i),Teff,R1,asemi,f,Ab

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
      subroutine getloggplanet(nunit,nplanet,npt,mp,np,mstar,rstar,seed)
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

c        mp(i)=pmass/(fourthirdsPi*prad*prad*prad)/1000.0d0
         mp(i)=log10(G*pmass/(prad*prad))+2.0d0 !add 2 to make cgs

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
c            goto 10 !we don't count negative mass
            pmass=(bc*ac*ac)**ct(1) !assume Mp << Ms
        endif

        prad=R1*rdr !planet radius
        
        mp(i)=pmass/(fourthirdsPi*prad*prad*prad)/1000.0d0
c         write(6,*) mp(i),rstar(k)*rdr*109.17,pmass/Mearth
c         read(5,*)
        
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
c            goto 10 !we don't count negative mass
            mp(i)=(bc*ac*ac)**ct(1)/Mearth !assume Mp << Ms
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
      
      read(10,*) dumr
      
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
      
