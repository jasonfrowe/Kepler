      program transitfit6
      implicit none
      integer nrad,ntheta,nfit,nunit,nplanet,i,nmax,npt,
     .  nptv,npta
      parameter(nfit=125,nmax=600000)
      integer dtype(nmax)
      real rbb(3,4)
      double precision sol(nfit),serr(nfit,2),err(nfit),time(nmax),
     .  flux(nmax),ferr(nmax),itime(nmax),Keplertime,vtime(nmax),
     .  vel(nmax),verr(nmax),vetime(nmax),kmag,kerr,aT(nmax),aM(nmax),
     .  aE(nmax),aIT(nmax),tmodel(nmax),bchi
C     Stuff for numerical recipes
C      integer ia(nfit)
c      double precision alpha(nfit,nfit),covar(nfit,nfit),dchistop
C     minpack
      integer iwa(nfit),lwa
      parameter(lwa=nmax*nfit+5*nmax*nfit)
      double precision fvec(nmax),wa(lwa),sol2(nfit)
      character*80 inputsol,obsfile,rvfile,cline
      character*3 titles(26)
c      common /Fitting/ nrad,ntheta,npta,nplanet,aIT
      common /Fitting2/ nrad,ntheta,nplanet,dtype,aT,aM,aE,aIT,sol,serr
      
      nrad=100   !number of rings to integrate star/planet surface 
      ntheta=100 !maximum number of angular divisions 
       
      if(iargc().lt.3) goto 901 !check number of commandline options
      
C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile) !read in arguement 
      nunit=10 !set file number
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      close(nunit)!release unit number as we are done with the file.

      kmag=0.0d0 !we can optionally read in Kepler-magnitude to set 
                 !photometric errors
      if(iargc().ge.4)then
        call getarg(4,cline)
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
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
      write(0,*) "done reading input solution"
c      call getfitpars(nunit,nfit,sol,serr,err)
      close(nunit) !release unit number as we are done with file
      
C     Store all the observations in the master file
      do 17 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
        aT(i)=time(i)
        aM(i)=flux(i)
        aE(i)=ferr(i)
        aIT(i)=itime(i)
c        write(0,*) "itime",itime(i)
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
      
      call pgopen('/xserve')
c     call pgopen('?')
      call pgpage()
      call PGPAP ( 15.0 ,1.0/3.0) 
      call pgsubp(3,1)

      call plotsetup(rbb,nfit,nplanet,sol,0)             
      call transitmodel(nrad,ntheta,nfit,nplanet,sol,npta,aT,aIT,dtype,
     .  tmodel,rbb,1)
      call plotphase(nfit,1,sol,npta,aT,aM,aE,dtype,tmodel,rbb)

C     Old crap ass NR routine that sucks balls..
c      call fittransitmodel(npta,aT,aM,aE,aIT,dtype,nfit,nplanet,sol,
c     .  serr,ia,covar,alpha,dchistop,err,rbb,tmodel,nrad,ntheta)

C     wonderful minpack routine that's faster than your mom.
      call fittransitmodel2(npta,nfit,sol,sol2,serr,fvec,iwa,lwa,wa)
      call plotsetup(rbb,nfit,nplanet,sol,0)  
      call transitmodel(nrad,ntheta,nfit,nplanet,sol,npta,aT,aIT,dtype,
     .  tmodel,rbb,1)
      call plotphase(nfit,1,sol,npta,aT,aM,aE,dtype,tmodel,rbb)      
      call exportfit(nfit,nplanet,sol,serr,err,titles)
      
      bchi=0.0d0
      do 13 i=1,npta
        bchi=bchi+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
c        write(0,*) i,bchi,tmodel(i)
 13   continue
      write(0,*) "bchi",bchi,npta-1
      bchi=dble(npta-1)/bchi
      write(0,*) "bchi_f:",bchi
     
      do 12 i=1,npt
        write(6,*) aT(i),aM(i),tmodel(i)
 12   continue
      
      
      call pgclos()

      goto 999
 901  write(0,*) "Usage: transitfit6 <photfile> <rvfile> <fitpars>"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 903  write(0,*) "Error opening ",obsfile
      goto 999
 904  write(0,*) "Error opening ",rvfile
      goto 999
 999  end      
      