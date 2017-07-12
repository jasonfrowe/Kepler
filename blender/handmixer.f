C     My evoling Blender Code
C
C     Todo:
C           dilution 
C
      program handmixer
      implicit none
      integer iargc,nfit,nunit,nmodelmax,nmodel,i,j,npt,nmax,npars,it,
     .  itmax,inl(3)
      parameter(nfit=18,nmodelmax=1000,nmax=600000)
      integer dtype(nmax),ipars(nfit),ia(nfit)
      double precision sol(nfit),serr(nfit,2),Dpvary(nfit),err(nfit),
     .  doe,toff,massrange(2),massi,Z,tage(nmodelmax),tTeff(nmodelmax),
     .  tlogL(nmodelmax),trad(nmodelmax),trho(nmodelmax),chio,
     .  tdrhodt(nmodelmax),dm,sol2(nfit),serr2(nfit,2),chisquared,
     .  time(nmax),mag(nmax),merr(nmax),itime(nmax),tmodel(nmax),
     .  MOSTtime,pvary(nfit),inclmin,dchistop,p(nfit),alpha(nfit,nfit),
     .  covar(nfit,nfit),rchi,minrad,nl1(76,11,8),nl2(76,11,8),
     .  nl3(76,11,8),nl4(76,11,8),G,Rsun,Msun,logg,Teff,gmag,
     .  Mg,M,Pi,tPi,Tsun,L,dist,dmag,gmagblend,gmaglimit,F2dF1,dil,bchi,
     .  gbright
      character*80 inputsol,obsfile,inputchar
      logical loop
      integer inda,indz,indm,nline
      parameter (nline=150)
      parameter (inda=3)
      parameter (indz=11)
      parameter (indm=35)
      integer a_one,z_one,m_one
      real*8 avalue(inda)
      real*8 zvalue(indz)
      real*8 xmass(indm)
      common /grid/avalue,zvalue,xmass
      common /single/a_one,z_one,m_one
      common /Fitting/ ipars,npt,itime,sol2,pvary,inclmin
      
c      open(unit=15,file="logfile.err")
      
      G=6.67259E-8
      Rsun=6.9599E10
      Msun=1.989E33
      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi      
      Mg=5.11 !absolute magnitude of the Sun in g-band
      Tsun=5781.6d0 !Solar temperature (K) 

C     Stuff to move to command line
      Teff=5500.0d0  !Temperture of KOI
c      gmag=11.864d0  !gmag of KOI
      gbright=0.0d0 !offset if we want to put a brighter star in
C     Define mass range
      massrange(1)=0.4 !range of Y^2 tracks
      massrange(2)=5.2
      dm=0.2  !mass interval to explore (not used)
      
      if(iargc().lt.5) goto 901 !check number of command line arguements

      call getarg(1,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
      close(nunit) !release unit number as we are done with file
      
C     Parse the name of the observations data file from the commandline
      call getarg(2,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      call readkeplc(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      close(nunit)!release unit number as we are done with the file.
      
      call getarg(3,inputchar)
      read(inputchar,*) massi
      if((massi.lt.massrange(1)).or.(massi.gt.massrange(2))) goto 901
      
      call getarg(4,inputchar)
      read(inputchar,*) gmag
      
      call getarg(5,inputchar)
      read(inputchar,*) gmagblend
      
      do 13 i=1,npt
        dtype(i)=0
 13   continue

      L=sol(3)*sol(3)*(Teff/Tsun)**4!luminosity of KOI
      M=Mg-2.5*log10(L)
c      write(0,*) "L: ",L,M
c      read(5,*)
      
C     read in limbdarkening table
      call limbdarkensetup(nl1,nl2,nl3,nl4)      
      
C     Get Chi-squared for 'best' solution
      chio=chisquared(npt,tmodel,time,mag,merr,itime,dtype,nfit,sol)
      chio=dble(npt-1)/chio !scale reduced chi-squared to 1 
      write(0,*) "BChi",chio
      
      itmax=1 !just in case...
      dchistop=1.0d-10 !criteria for good fit.
c      massi=massrange(1) !now defined from commandline 
      Z=0.02
      gmaglimit=25.0d0
      dmag=1.0 !how much fainter to make additional star
      
C     We will now run for one mass at a time.
c      do 11 while(massi.lt.massrange(2))
      
c        gmagblend=gmag-gbright !reset gmagblend
        
C       Get evolution track       
        call gettrack(massi,Z,nmodelmax,nmodel,tage,tTeff,tlogL,
     .      trad,trho,tdrhodt)
        
c        do 15 while(gmagblend.lt.gmaglimit)
            
            F2dF1=10**((gmagblend-gmag)/-2.5d0)
            dil=1.0d0/(1.0d0+F2dF1)
c            read(5,*)
c            goto 16
     
C       get minimum radius
C       Get limbdarkening
            if(tTeff(1).le.13000.0)then
                inl(1)=(int(tTeff(1))-3500)/250+1
            elseif(tTeff(1).gt.13000)then
                inl(1)=(int(tTeff(1))-13000)/1000+39
            endif
            logg=log10(G*massi*Msun/(Rsun*Rsun*trad(1)*trad(1)))
            inl(2)=int(10.0*logg)/5+1      
            inl(3)=6 ![Fe/H]=0

C       Copy original solution
            do 14 j=1,nfit
                sol2(j)=sol(j)
                serr2(j,1)=0.0d0
                serr2(j,2)=0.0d0
 14         continue
 
            serr2(2,2)=0.0 !fit companion mass
            serr2(9,2)=0.0 !fit albedo
            serr2(16,2)=0.0 !fit occcultation
 
            serr2(4,2)=-1.0 !fit companion radius
            serr2(5,2)=-1.0 !fit period
            serr2(3,2)=-1.0 !fit inclination
            serr2(7,2)=-1.0 !fit T0
            serr2(8,2)=-1.0 !fit zpt
            sol2(6)=90.00
            sol2(1)=massi
            sol2(3)=trad(1)
            sol2(4)=sol(4)*sqrt(1.0d0/(1.0d0-dil))
C       update limb-darkening
            sol2(10)=nl1(inl(1),inl(2),inl(3))
            sol2(11)=nl2(inl(1),inl(2),inl(3))
            sol2(12)=nl3(inl(1),inl(2),inl(3))
            sol2(13)=nl4(inl(1),inl(2),inl(3))
            sol2(18)=dil
            
c            write(15,500) massi,gmagblend,dil,sol2(3),sol2(4)/10.05
            if(sol2(4)/10.05d0.ge.sol2(3)) goto 17
            
            it=0
            loop=.true.
            do while(loop) 
                it=it+1
                call fittransitmodel2(npt,time,mag,merr,dtype,nfit,
     .              sol2,serr2,npars,ipars,p,pvary,Dpvary,dchistop,
     .              alpha,covar,ia,inclmin,err)
                if(it.gt.itmax)loop=.false.
            enddo
            minrad=sol2(3)
 17         continue !jump here if companion radius is larger.

            do 10 i=1,nmodel
                if((tage(i).gt.16.0d0).or.(trad(i).lt.minrad)) goto 10
            
C           Get limbdarkening
                if(tTeff(i).le.13000.0)then
                    inl(1)=(int(tTeff(i))-3500)/250+1
                elseif(tTeff(i).gt.13000)then
                    inl(1)=(int(tTeff(i))-13000)/1000+39
                endif
                logg=log10(G*massi*Msun/(Rsun*Rsun*trad(i)*trad(i)))
                inl(2)=int(10.0*logg)/5+1      
                inl(3)=6 ![Fe/H]=0
            
C           Copy original solution
                do 12 j=1,nfit
                    sol2(j)=sol(j)
                    serr2(j,1)=0.0d0
                    serr2(j,2)=0.0d0
 12             continue
 
                serr2(2,2)=0.0 !fit companion mass
                serr2(9,2)=0.0 !fit albedo
                serr2(16,2)=0.0 !fit occcultation
                
                serr2(4,2)=-1.0 !fit companion radius
                serr2(5,2)=-1.0 !fit period
                serr2(6,2)=-1.0 !fit inclination
                serr2(7,2)=-1.0 !fit T0
                serr2(8,2)=-1.0 !fit zpt
                sol2(6)=89.85
                sol2(1)=massi
                sol2(3)=trad(i)
                sol2(4)=sol(4)*sqrt(1.0d0/(1.0d0-dil))
C           update limb-darkening
                sol2(10)=nl1(inl(1),inl(2),inl(3))
                sol2(11)=nl2(inl(1),inl(2),inl(3))
                sol2(12)=nl3(inl(1),inl(2),inl(3))
                sol2(13)=nl4(inl(1),inl(2),inl(3))
                sol2(18)=dil
c            sol(18)=10**tlogL(i)
c            write(0,*) massi,trad(i),tlogL(i)
 
                if(sol2(4)/10.05d0.ge.sol2(3)) goto 18

                it=0
                loop=.true.
                do while(loop) 
                    it=it+1
                    call fittransitmodel2(npt,time,mag,merr,dtype,nfit,
     .                  sol2,serr2,npars,ipars,p,pvary,Dpvary,dchistop,
     .                  alpha,covar,ia,inclmin,err)
                    if(it.gt.itmax)loop=.false.
                enddo
 18             continue !jump here if companion mass greater
                rchi=chisquared(npt,tmodel,time,mag,merr,
     .              itime,dtype,nfit,sol2)*chio
                write(6,500) massi,trad(i),sol2(4),sol2(6),tTeff(i),
     .              tlogL(i),gmagblend,dil,rchi
 10         continue
c 16         gmagblend=gmagblend+dmag !iterate to fainter star
c 15     continue
c        massi=massi+dm !move to next mass track
c 11   continue
 500  format(10(1X,1PE17.10)) 
      
      goto 999
 901  write(0,*) "Usage: handmixer f1.dat kld.d.dat smass gmag gmagblend
     ."
      write(0,*) "smass much be between 0.4 and 5.2"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 903  write(0,*) "Cannot open ",obsfile
      goto 999
 999  continue
c      close(15) !close error log file
      end