C     This program allows one to explore parameter space to find a good
C     first guess for fitting a transit
      program modelbuild
      implicit none
C     nmax: the maximum number of data points that can be read in
      integer nmax,iargc,nunit,npt,nfit,nexpl,nptv,i,npta
      parameter(nmax=650000,nfit=18)
      integer dtype(nmax)
      real px(nmax),py(nmax),pz(nmax)
      double precision time(nmax),mag(nmax),merr(nmax),itime(nmax),
     .  MOSTtime,bb(4),bbp(4),bbp2(4),sol(nfit),serr(nfit,2),
     .  Dpvary(nfit),toff,phase(nmax),tbin,avgitime,mean,tmodel(nmax),
     .  sol2(nfit),err(nfit),doe,vtime(nmax),vel(nmax),verr(nmax),
     .  vetime(nmax),bbp3(4),avgvtime,aT(nmax),aM(nmax),aE(nmax),
     .  aIT(nmax)
      character*3 titles(nfit)
      character*80 inputsol,obsfile,command,rvfile
      logical loop
       
      loop=.true.  !initialze loop 
      nexpl=0      !if nexpl.eq.0, then we are going to explore chi-sq
      
      if(iargc().lt.1) goto 901 !check number of command line arguements
      call getarg(1,obsfile) !get filename for input solution
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      call readkeplc(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      do 10 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
        aT(i)=time(i)
        aM(i)=mag(i)
c        merr(i)=0.0003
        aE(i)=merr(i)
        aIT(i)=itime(i)
 10   continue
      npta=npt
      close(nunit)!release unit number as we are done with the file.

C     We can bin the data
      tbin=30.0
c      call bindt(npt,time,mag,merr,itime,tbin)

      if(iargc().ge.2) then
        call getarg(2,rvfile) !get filename for RV data
        if(rvfile.eq.'null')then
            nptv=0
        else
            nunit=10
            open(unit=nunit,file=rvfile,status='old',err=904)
            call readrv(nunit,nmax,nptv,vtime,vel,verr,vetime)
            close(nunit)
c            write(6,*) "rv:",(vtime(i),i=1,nptv)
        endif
        do 11 i=1,nptv
            dtype(npt+i)=1 !mark data as RV
            aT(npt+i)=vtime(i)
            aM(npt+i)=vel(i)
            aE(npt+i)=verr(i)
            aIT(npt+i)=vetime(i)
 11     continue
        npta=npta+nptv
      endif
      
      if(iargc().ge.3) then
        call getarg(3,inputsol) !get filename for input solution
        nunit=10 !unit number used for file input
        open(unit=nunit,file=inputsol,status='old',err=902)
C       We start by reading in solution from input file
        call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
        close(nunit) !release unit number as we are done with file
      else
C     Import a default solution set to start with
        call defaultpars(nfit,sol,serr,Dpvary,err,doe,toff)
      endif
      
      call opengraphics()
      call pgslw(2)
      
      do while(loop)
        call pgeras()
C       Plot of the observations
        call plotdata(npt,time,mag,merr,-1.0d0,-1.0d0,1,MOSTtime,nmax,
     .      px,py,pz,bb)
c       Phase plot of the observations
        call plotph(npt,time,mag,merr,sol(5),0.0d0,1.0d0,2,-1,nmax,px,
     .      py,pz,phase,toff,1,bbp)
     
        call plotph(npt,time,mag,merr,sol(5),0.700d0,0.800d0,3,-1,nmax,
     .      px,py,pz,phase,toff,1,bbp2)
     
        if(nptv.gt.0)then
        call plotrvph(nptv,vtime,vel,verr,sol(5),0.0d0,1.0d0,4,17,nmax,
     .      px,py,pz,phase,toff,0,bbp3)
        endif

C       Plotting the solution     
        avgitime=mean(npt,itime)
        avgvtime=mean(nptv,vetime)
        call modelplot(nfit,sol,nmax,phase,bb,bbp,bbp2,bbp3,avgitime,
     .      avgvtime,toff,nptv)
     
        call printsol(nfit,sol,serr,Dpvary,err,toff)

        write(0,*) "Enter command: "
        read(5,500) command
 500    format(A80)
        call parsecommand(command,nfit,sol,toff,loop,nexpl)
        if(nexpl.ne.0) then
            call explore(nexpl,npta,aT,aM,aE,aIT,dtype,tmodel,nfit,sol,
     .          sol2,serr,Dpvary)
            nexpl=0
        endif
      enddo
           
      call pgslw(1)
      call closegraphics()
      
      call exportfit(nfit,sol,serr,Dpvary,err,doe,toff,titles)
      
      goto 999
 901  write(0,*) "Usage: modelbuild <obsfile> <rvfile> <modelpars>"
      write(0,*) " <obsfile>: file containing photometry"
      write(0,*) " <rvfile> : file containing RV photometry"
      write(0,*) "            enter 'null' to proceed without RV"
      write(0,*) " <modelpars>: model parameters to adjust"
      write(0,*) "              omit to use defaults"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 903  write(0,*) "Cannot open ",obsfile
      goto 999
 904  write(0,*) "Cannot open ",rvfile
      goto 999
 999  end




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine parsecommand(command,nfit,sol,toff,loop,nexpl)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nexpl
      double precision sol(nfit),toff
      character*80 command
      logical loop
      
      if(command(1:3).eq."sma") then
        read(command(5:80),*) sol(1)
      elseif(command(1:3).eq."pma") then
        read(command(5:80),*) sol(2)
      elseif(command(1:3).eq."sra") then
        read(command(5:80),*) sol(3)
      elseif(command(1:3).eq."pra") then
        read(command(5:80),*) sol(4)
      elseif(command(1:3).eq."per") then
        read(command(5:80),*) sol(5)
      elseif(command(1:3).eq."inc") then
        read(command(5:80),*) sol(6)
      elseif(command(1:3).eq."epo") then
        read(command(5:80),*) sol(7)
      elseif(command(1:3).eq."zpt") then
        read(command(5:80),*) sol(8)
      elseif(command(1:3).eq."alb") then
        read(command(5:80),*) sol(9)
      elseif(command(1:3).eq."nl1") then
        read(command(5:80),*) sol(10)
      elseif(command(1:3).eq."nl2") then
        read(command(5:80),*) sol(11)
      elseif(command(1:3).eq."nl3") then
        read(command(5:80),*) sol(12)
      elseif(command(1:3).eq."nl4") then
        read(command(5:80),*) sol(13)
      elseif(command(1:3).eq."nla") then
        read(command(5:80),*) sol(10),sol(11),sol(12),sol(13)
      elseif(command(1:3).eq."ecw") then
        read(command(5:80),*) sol(14)
      elseif(command(1:3).eq."esw") then
        read(command(5:80),*) sol(15)
      elseif(command(1:3).eq."ted") then
        read(command(5:80),*) sol(16)
      elseif(command(1:3).eq."vof") then
        read(command(5:80),*) sol(17)
      elseif(command(1:3).eq."dil") then
        read(command(5:80),*) sol(18)
      elseif(command(1:3).eq."off") then
        read(command(5:80),*) toff
      elseif(command(1:3).eq."end") then
        loop=.false.
      elseif(command(1:3).eq."exp") then
        read(command(5:80),*) nexpl
      endif
        
      
      return
      end
      
