CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program transitmcmc5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Fits for stellar density (a/R*)
C     Jason Rowe - jasonfrowe@gmail.com
      implicit none
      integer iargc,nfit,npt,nmax,i,nunit,nptv,npta,j,naccept,ng,
     .  nplanet,nfitm,nbuffer,nmov,nup,nb,npars,nupdate,nplanetmax
      parameter(nfitm=108,nmax=2000000,nbuffer=500,nplanetmax=10)
      integer dtype(nmax),niter,flag,nacor,
     .  nacorsub,ngcor(nfitm),ngcorsub(nfitm),
     .  naprob,naprobsub,ngprob(nfitm),
     .  ngprobsub(nfitm),ntt(nplanetmax),ngs(nfitm),nas
      double precision sol(nfitm),time(nmax),dt,tmodel(nmax),
     .  flux(nmax),ferr(nmax),exptime(nmax),Keplertime,serr(nfitm,2),
     .  err(nfitm),vtime(nmax),vel(nmax),verr(nmax),vetime(nmax),
     .  aT(nmax),aM(nmax),aE(nmax),aIT(nmax),kmag,kerr,ran2,dumr,
     .  sol2(nfitm),rchi,bchi,gasdev,accrate,dil(2),
     .  buffer(nfitm,nbuffer),corscale,gscale(nfitm),gratio(nfitm),
     .  fchi,chiold
C     Random number generation
      integer now(3),seed
      character*80 inputsol,obsfile,rvfile,cline,chseed,rhofile,ttfile
C     TT variations
      double precision tobs(nplanetmax,nmax),omc(nplanetmax,nmax)      
      integer nfrho
      double precision rhoread(8),rhoi,rhoierr(9),rhoin(9)
c      data rhoread /1143.834961,  220.177612,  169.443970, -280.480774, 
c     .  339.912354, -545.593018,  550.404907, -717.149841/
      data rhoin/-1.0d2,-3.0d0,-2.0d0,-1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,
     .  1.0d2/      
      
      do 23 i=1,nfitm
        gscale(i)=1.0d0
        ngcor(i)=0
        ngcorsub(i)=0
        ngprob(i)=0
        ngprobsub(i)=0
 23   continue
      
      corscale=1.0d0
      nacor=0
      nacorsub=0
      naprob=0
      naprobsub=0
     
      nmov=0 !initialize buffer size for dmcmc
      nup=0 !when to update buffer (0=no,1=yes)
      nb=0  !which buffer element to replace (oldest)
      nupdate=0 !number of times scales have been updated
      
      if(iargc().lt.4) goto 901

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,flux,ferr,exptime,Keplertime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,exptime,Keplertime)
      close(nunit)!release unit number as we are done with the file.
C     apply a boxcar filter
c      boxbin=1.0 !filter width (days)
c      call boxcar(npt,time,mag,merr,boxbin)

C     Get density constraints..
      nfrho=1
      if(iargc().ge.5)then
        call getarg(5,rhofile)
        if(rhofile.eq.'null') goto 25
        open(unit=nunit,file=rhofile,status='old',err=905)
        read(nunit,*) (rhoread(i),i=1,8)
        close(nunit)
        do 26 i=1,8
            rhoread(i)=rhoread(i)*1000.0d0 !convert to kg/m^3
 26     continue
        rhoi=rhoread(1)
        rhoierr(1)=-rhoi
        rhoierr(2)=rhoread(8)
        rhoierr(3)=rhoread(6)
        rhoierr(4)=rhoread(4)
        rhoierr(5)=0.0d0
        rhoierr(6)=rhoread(3)
        rhoierr(7)=rhoread(5)
        rhoierr(8)=rhoread(7)
        rhoierr(9)=rhoierr(8)*10.0d0
        if(rhoread(2).gt.0.0d0) then 
            nfrho=0
            write(0,*) rhoi,(rhoierr(i),i=1,9)
            write(0,*) "density constraints are ON"
        endif
 25     continue
      endif
      if(nfrho.eq.1) write(0,*) "density constraints are OFF"

      kmag=0.0d0
      if(iargc().ge.6)then
        call getarg(6,cline)
        read(cline,*) kmag
      endif
      if(kmag.gt.0.0)then !estimate of expected scatter
C       quadratic version
c        kerr=(4.0d0*(kmag-8.0d0)**2.0+10.0d0)*1.0d-6
C       Powerlaw version
         kerr=(10.0d0**(0.23*kmag-0.9))*1.0d-6
         do 11 i=1,npt
            ferr(i)=kerr
 11      continue
      else
        kerr=100.0*1.0d-6
      endif 
      write(0,*) "KMAG,KERR:",kmag,kerr
      


      call getarg(2,rvfile) !get filename for RV data
      if(rvfile.eq.'null')then
        nptv=0
      else
        nunit=10
        open(unit=nunit,file=rvfile,status='old',err=904)
        call readrv(nunit,nmax,nptv,vtime,vel,verr,vetime)
        close(nunit)
c        write(6,*) "rv:",(vtime(i),i=1,nptv)
      endif

      call getarg(3,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfitm,nplanet,sol,serr,err)
      close(nunit) !release unit number as we are done with file
      write(0,*) "nPlanet: ",nplanet
      nfit=nplanet*10+8
      
c      write(0,*) (serr(i,2),i=1,nfit)
    
C     Store all the observations in the master file
      do 17 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
        aT(i)=time(i)
        aM(i)=flux(i)
        aE(i)=ferr(i)
        aIT(i)=exptime(i)
 17   continue
      npta=npt
      do 18 i=1,nptv
        dtype(npt+i)=1 !mark data as RV
        aT(npt+i)=vtime(i)
        aM(npt+i)=vel(i)
        aE(npt+i)=verr(i)
        aIT(npt+i)=vetime(i)
 18   continue
      npta=npta+nptv

      do 24 i=1,nplanet
        if(iargc().ge.6+i)then
            call getarg(6+i,ttfile)

            if(ttfile.eq.'null')then
                ntt(i)=0
            else
                nunit=10
                open(unit=nunit,file=ttfile,status='old',err=905)
                call readttfile(nunit,nplanetmax,nmax,i,ntt,tobs,omc)
                close(nunit)
                write(0,*) "Using TT files"
            endif
            
        else
            ntt(i)=0
        endif
c        write(0,*) "ntt",i,ntt(i)
 24   continue

      
      write(0,*) "Initialization of random number"
C     Random number initialization(so we get a different seed each time)
c      if(iargc().ge.7)then
c        call getarg(7,chseed)
c        read(chseed,*) seed
c      else
        call itime(now)
        seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
c      endif
      write(0,*) "Seed: ",seed
      dumr=ran2(-seed)
      
      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npta,aT,aIT,
     .  ntt,tobs,omc,tmodel,dtype)
      
      bchi=0.0d0
      do 12 i=1,npta
        bchi=bchi+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
 12   continue
      fchi=bchi !keep tabs on the best chi-squared
      write(0,*) "bchi",bchi,npta-1
      bchi=dble(npta-1)/bchi
      write(0,*) "bchi_f:",bchi
c      bchi=1.0d0

C     Save dilution parameters
      dil(1)=sol(6)
      dil(2)=serr(6,2)
c      write(0,*) "dil",dil(1),dil(2)

      nas=0 !initalized flag for updating jump scale
      npars=1 !counting number parameters that are fitted
C     Initialization of Markov chain
      do 13 i=1,nfit
         ngs(i)=0 ! initialize flag for updating jump scale
c        sol(i)=serr(i,1)+serr(i,2)*(2.0d0*ran2(seed)-1.0d0)

C     Eccentricity constraints
c        if(serr(i,1).eq.0.0d0) then
c            err(i,1)=sol(i)
c            err(i,2)=sol(i)
c        endif

        if((i.lt.3).or.(i.gt.6))then !skip dilution and do limb-darkening together
c            if(serr(i,2).gt.0.0) npars=npars+1 !counting fitted pars
 22         sol(i)=gasdev(seed)*serr(i,2)+sol(i)

            if(i.eq.2) then !limb-darkening updates.
               sol(3)=gasdev(seed)*serr(3,2)+sol(3)
               sol(4)=gasdev(seed)*serr(4,2)+sol(4)
               sol(5)=gasdev(seed)*serr(5,2)+sol(5)
               if((sol(2).eq.0.0).and.(sol(3).eq.0.0))then  !Kipping Limb-darkening
                  if((sol(4).lt.0.0).or.(sol(4).gt.1.0).or.
     .             (sol(5).lt.0.0).or.(sol(5).gt.1.0))then
                     goto 22  !if invalid, draw again
                  endif
               elseif((sol(4).eq.0.0).and.(sol(5).eq.0.0))then !Quadratic limb-darkening
                  if((sol(2)+sol(3).gt.1.0).or.(sol(2).lt.0).or.
     .             (sol(2)+2.0d0*sol(3).lt.0))then
                     goto 22 !if invalid, draw again
                  endif
               endif
            endif

            if(i.gt.9)then
                j=i-10*((i-9)/10)
                if((j.eq.11).and.(sol(i).lt.0.0d0)) goto 22
C     This line keeps b < 1 (at least to start)
c                if((j.eq.11).and.(sol(i).gt.1.0d0)) goto 22
                if((j.eq.12).and.(sol(i).lt.0.0d0)) goto 22
c     This line keeps K > 0
c                if((j.eq.15).and.(sol(i).lt.0.0d0)) goto 22
            endif

        endif
 13   continue
      sol(1)=abs(sol(1)) !mean-stellar density must be positive!
 21   sol(6)=gasdev(seed)*dil(2)+dil(1)
      if(sol(6).lt.0.0d0) goto 21

      chiold=0.0 !keeps track of previous chi-sq.
      write(6,*) nfit
      naccept=0 !calculate accept rate
      call getarg(4,cline)
      read(cline,*) niter
      do 19 j=1,niter
        call dmcmc(npta,aT,aM,aE,aIT,dtype,nfitm,nfit,nplanet,sol,sol2,
     .      serr,err,tmodel,rchi,seed,flag,bchi,ng,dil,nup,nb,nmov,
     .      nbuffer,buffer,npars,corscale,nacor,nacorsub,naprob,
     .      naprobsub,gscale,ngcor,ngcorsub,ngprob,ngprobsub,nupdate,
     .      gratio,nfrho,rhoi,rhoierr,rhoin,nplanetmax,nmax,ntt,tobs,
     .      omc,ngs,nas,chiold)
        if(flag.eq.0)then
            do 20 i=1,nfit
                sol(i)=sol2(i)

C              Eccentricity constraints
c                if(serr(i,1).eq.0.0d0) then
c                  err(i,1)=sol(i)
c                  err(i,2)=sol(i)
c                endif

 20         continue
            naccept=naccept+1
        endif
        if(j.gt.1) write(6,500) rchi/bchi,flag,ng,
     .      (sol(i),i=1,nplanet*10+8)
 500    format(1PE17.10,2(1X,I2),108(1X,1PE17.10))

c        write(0,*) rchi,rchi/bchi,fchi
        if(rchi/bchi.lt.fchi)then  !check for new chi-squared min
            fchi=rchi/bchi
            bchi=dble(npta-1)/fchi
            write(0,*) "fchi",fchi,bchi
        endif
  
 19   continue
      accrate=dble(naccept)/dble(niter)
      write(0,*) "Accept Rate: ",accrate
      
      goto 999
 901  write(0,*) "Usage: transitfit5 <photfile> <rvfile> <fitpars> <nite
     .r> <rhostar> <kmag>"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 903  write(0,*) "Error opening ",obsfile
      goto 999
 904  write(0,*) "Error opening ",rvfile
      goto 999
 905  write(0,*) "Error opening ",rhofile
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ovrwrt (line, iwhich)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Taken from the DAOPhot code by Stetson.
C     This is a cool routine for writing lots of info to the screen
C     and keeping everything on the same line.
C
      character*(*) line
      character*79 output
      integer len
      if (iwhich .eq. 1) then
         write (0,1) line
    1    format (a)
      else if (iwhich .eq. 2) then
         if (len(line) .lt. 79) then
            output = ' '
            output = line
            write (0,2) output, char(13), char(13)
            write (0,2) output, char(13), char(13)
            write (0,2) output, char(13), char(13)
    2       format (a, 2a1, $)
         else
            write (0,2) line, char(13), char(13)
         end if
      else if (iwhich .eq. 3) then
         write (0,3) line
    3    format (a)
      else
         write (0,4) line, char(13), char(13)
    4    format (/a, 2a1, $)
         write (0,2) line, char(13), char(13)
      end if
      return
      end
