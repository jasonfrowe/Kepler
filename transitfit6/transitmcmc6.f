      program transitfit6
      implicit none
      integer nrad,ntheta,nfit,nunit,nplanet,i,nmax,npt,nptv,npta,nfitm,
     .  j
      parameter(nfitm=125,nmax=600000)
      integer dtype(nmax)
      real rbb(3,4)
      double precision sol(nfitm),serr(nfitm,2),err(nfitm),time(nmax),
     .  flux(nmax),ferr(nmax),exptime(nmax),Keplertime,vtime(nmax),
     .  vel(nmax),verr(nmax),vetime(nmax),kmag,kerr,aT(nmax),aM(nmax),
     .  aE(nmax),aIT(nmax),tmodel(nmax)
      character*80 inputsol,obsfile,rvfile,cline
      character*3 titles(26)
C     mcmc vars
      integer nfrho,ngcor(nfitm),ngcorsub(nfitm),ngprob(nfitm),
     .  ngprobsub(nfitm),nacor,nacorsub,naprob,naprobsub,seed,now(3),
     .  npars,naccept,niter,flag,ng,nup,nb,nmov,nbuffer,nupdate
      parameter (nbuffer=500)
      double precision rhoread(8),rhoi,rhoierr(9),rhoin(9),gasdev,
     .  gscale(nfitm),corscale,dumr,ran2,bchi,fchi,dil(2),sol2(nfitm),
     .  rchi,buffer(nfitm,nbuffer),gratio,accrate
      character*80 rhofile,chseed
      data rhoin/-1.0d2,-3.0d0,-2.0d0,-1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,
     .  1.0d2/      

      nrad=100   !number of rings to integrate star/planet surface 
      ntheta=100 !maximum number of angular divisions 

C     Initialize variables that control adoptive Gibbs sampling
      do 23 i=1,nfitm
        gscale(i)=1.0d0 !Gibbs scale factor
        ngcor(i)=0
        ngcorsub(i)=0
        ngprob(i)=0
        ngprobsub(i)=0
 23   continue
      corscale=1.0d0 !vector scale factor
      nacor=0
      nacorsub=0
      naprob=0
      naprobsub=0
      nmov=0 !initialize buffer size for dmcmc
      nup=0 !when to update buffer (0=no,1=yes)
      nb=0  !which buffer element to replace (oldest)
      nupdate=0 !number of times scales have been updated
       
      if(iargc().lt.4) goto 901 !check number of commandline options
      
C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile) !read in arguement 
      nunit=10 !set file number
      open(unit=nunit,file=obsfile,status='old',err=903)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,exptime,Keplertime)
      close(nunit)!release unit number as we are done with the file.

C     Get density constraints..
      nfrho=1
      if(iargc().ge.5)then
        call getarg(5,rhofile)
        if(rhofile.eq.'null') goto 25
        open(unit=nunit,file=rhofile,status='old',err=905)
        read(nunit,*) (rhoread(i),i=1,8)
        close(nunit)
        do 26 i=1,8
            rhoread(i)=rhoread(i)*1000.0d0
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
            write(0,*) "density constraints are ON"
        endif
 25     continue
      endif
      if(nfrho.eq.1) write(0,*) "density constraints are OFF"

      kmag=0.0d0 !we can optionally read in Kepler-magnitude to set 
                 !photometric errors
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
 11     continue
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
      write(0,*) "reading in input solution"
      call getfitpars(nunit,nfitm,nplanet,sol,serr,err)
      write(0,*) "done reading input solution"
c      call getfitpars(nunit,nfit,sol,serr,err)
      close(nunit) !release unit number as we are done with file
      write(0,*) "nPlanet: ",nplanet
      nfit=nplanet*11+15
      
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

      write(0,*) "Initialization of random number"
C     Random number initialization(so we get a different seed each time)
      if(iargc().ge.7)then
        call getarg(7,chseed)
        read(chseed,*) seed
      else
        call itime(now)
        seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      endif
      write(0,*) "Seed: ",seed
      dumr=ran2(-seed)
      
      call pgopen('/xserve')
c     call pgopen('?')
      call pgpage()
      call PGPAP ( 15.0 ,1.0/3.0) 
      call pgsubp(3,1)

      call plotsetup(rbb,nfit,nplanet,sol,0)             
      call transitmodel(nrad,ntheta,nfitm,nplanet,sol,npta,aT,aIT,dtype,
     .  tmodel,rbb,1)
      call plotphase(nfit,1,sol,npta,aT,aM,aE,dtype,tmodel,rbb)

      bchi=0.0d0
      do 12 i=1,npta
        bchi=bchi+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
c        write(0,*) i,bchi,tmodel(i)
 12   continue
      fchi=bchi !keep tabs on the best chi-squared
      write(0,*) "bchi",bchi,npta-1
      bchi=dble(npta-1)/bchi
      write(0,*) "bchi_f:",bchi

C     Save dilution parameters
      dil(1)=serr(6,1)
      dil(2)=serr(6,2)

      npars=1 !controls how much history to keep
C     Initialization of Markov chain
c      write(0,*) "nfit:",nfit
      do 13 i=1,nfit
        if(i.ne.6)then
 22         sol2(i)=gasdev(seed)*serr(i,2)+sol(i)
            if((i.eq.1).and.(sol2(i).lt.0.0)) goto 22
            if((i.eq.12).and.(sol2(i).lt.0.0)) goto 22
            if((i.eq.12).and.(sol2(i).ge.1.0)) goto 22          
            if(i.gt.15)then
                j=i-11*((i-15)/11)
                if((j.eq.18).and.(sol2(i).lt.0.0d0)) goto 22
                if((j.eq.19).and.(sol2(i).lt.0.0d0)) goto 22
                if((j.eq.22).and.(sol2(i).lt.0.0d0)) goto 22
            endif
        endif
c        write(0,*) i,sol2(i)
 13   continue
      do 14 i=1,nfit
        sol(i)=sol2(i)
 14   continue
 21   sol(6)=gasdev(seed)*dil(2)+dil(1)
      if(sol(6).lt.0.0d0) goto 21

      write(6,*) nfit
      rchi=-1.0d0 !initialize chi-square for MCMC
      naccept=0 !calculate accept rate
      call getarg(4,cline)
      read(cline,*) niter
      do 19 j=1,niter
        call dmcmc(npta,aT,aM,aE,aIT,dtype,nfitm,nfit,nplanet,sol,sol2,
     .      serr,err,tmodel,rchi,seed,flag,bchi,ng,dil,nup,nb,nmov,
     .      nbuffer,buffer,npars,corscale,nacor,nacorsub,naprob,
     .      naprobsub,gscale,ngcor,ngcorsub,ngprob,ngprobsub,nupdate,
     .      gratio,nfrho,rhoi,rhoierr,rhoin,nrad,ntheta,rbb)
        if(flag.eq.0)then
            do 20 i=1,nfit
                sol(i)=sol2(i)
 20         continue
            naccept=naccept+1
        endif
        if(j.gt.1) write(6,500) rchi/bchi,flag,ng,
     .      (sol(i),i=1,nplanet*11+15)
 500    format(1PE17.10,2(1X,I2),125(1X,1PE17.10))

c        write(0,*) rchi,rchi/bchi,fchi
        if(rchi/bchi.lt.fchi)then  !check for new chi-squared min
            fchi=rchi/bchi
            bchi=dble(npta-1)/fchi
            write(0,*) "fchi",fchi,bchi
        endif
  
 19   continue
      accrate=dble(naccept)/dble(niter)
      write(0,*) "Accept Rate: ",accrate

c      call transitmodel(nrad,ntheta,nfit,nplanet,sol,npta,aT,aIT,dtype,
c     .  tmodel,rbb,1)
c      call plotphase(nfit,1,sol,npta,aT,aM,aE,dtype,tmodel,rbb)      
c      call exportfit(nfit,nplanet,sol,serr,err,titles)      
      
      call pgclos()

      goto 999
 901  write(0,*) "Usage: transitmcmc6 <photfile> <rvfile> <fitpars> <nit
     .er> [rhoboot.dat] [kmag] [seed]"
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
      